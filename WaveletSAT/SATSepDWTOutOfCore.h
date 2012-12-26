#pragma once

#include <map>	

#include <vector>
using namespace std;
#include <math.h>

#include "SepDWTHeader.h"
#include "SATSepDWTNetCDF.h"
#include "SepDWTPool.h"

#include "liblog.h"	
#include "libbuf.h"

#if	WITH_NETCDF
#include <netcdf.h>
#include "lognc.h"
#endif	// #if	WITH_NETCDF

namespace WaveletSAT
{
	template<typename T>
	// The class that access the coefficients from files (in NetCDF format)
	class CSATSepDWTOutOfCore:
		public CSATSepDWTNetCDF,
		public CSepDWTHeader
	{
protected:	
			const char* szFilepath;	

			vector<size_t>			vuHeaderOffset;
			vector<unsigned short>	vusHeaderCount;

			//! #Coefs stored in full arrays
			size_t uNrOfCoefsInFullArray;	

			//! vector of pointers to the coefficient pools. 
			/*!
			If the poiner is NULL, it means that the coefficients are out of core.
			*/
			vector< CSepDWTPool<T, unsigned short>* > vpcCoefPools;

			TBuffer<TYPE_COEF_BIN>		pCoefBin;
			TBuffer<TYPE_COEF_VALUE>	pCoefValue;
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0a00,
			PARAMETER_END
		};

		virtual	
		void
		_SetLong(
			int eName,
			long lValue,
			void* _Reserved = NULL
		)
		{
			CSATSepDWTNetCDF::_SetLong(eName, lValue);
			CSepDWTHeader::_SetLong(eName, lValue);
		}

		virtual
		void
		_Allocate(
					void *_Reserved = NULL
					)
		{
			vuHeaderOffset.resize(uNrOfCoefs);
			vusHeaderCount.resize(uNrOfCoefs);
			pCoefBin.alloc(UGetNrOfBins());
			pCoefValue.alloc(UGetNrOfBins());

			// decide which basis should be in core
			uNrOfCoefsInFullArray = min(
				(size_t)floor( (double)uSizeOfFullArrays/(double)(UGetNrOfBins() * sizeof(T))), 
				uNrOfCoefs);

			/////////////////////////////////////////////////////////////////
			size_t uNrOfDims = this->UGetNrOfDims();
			size_t uLevelsProduct = 1;	// product of levels from all dim
			for(size_t d = 0; d < uNrOfDims; d++)
			{
				size_t uNrOfDimLevels = (size_t)ceil((log((double)vuDimLengths[d]) / M_LN2)) + 1;
				uLevelsProduct *= vuDimLevels[d];
			}
			vector<size_t> vuOptimalDimLevel;
			vuOptimalDimLevel.resize(uNrOfDims);
			
			if( uNrOfCoefsInFullArray )
			{
				size_t uDiff = uSizeOfFullArrays;
				double dCurrentAspectRatio = -1.0; // ADD-BY-LEETEN 10/31/2012
				for(size_t l = 0; l < uLevelsProduct; l++)
				{
					vector<size_t> vuLevel;
					_ConvertIndexToSub(l, vuLevel, vuDimLevels);
					size_t uSize = 1;

					double dMaxLength = -HUGE_VAL;
					double dMinLength = +HUGE_VAL;

					for(size_t d = 0; d < uNrOfDims; d++)
					{
					  size_t uLength = (size_t)(1 << vuLevel[d]);
					  uSize *= uLength;
					  dMaxLength = max(dMaxLength, (double)uLength);
					  dMinLength = min(dMinLength, (double)uLength);
					}
					double dAspectRatio = dMaxLength / dMinLength;

					if( uSize <= uNrOfCoefsInFullArray )
					{
						size_t uNewDiff = uNrOfCoefsInFullArray - uSize;
						if( uNewDiff <= uDiff )
						  if( uNewDiff < uDiff ||
						      (dCurrentAspectRatio < 0.0 || dAspectRatio <= dCurrentAspectRatio ) )
						    {
							vuOptimalDimLevel = vuLevel;
							uDiff = uNewDiff;
							dCurrentAspectRatio = dAspectRatio;
						    }
					}
				}
			}

			size_t uNrOfOutOfCores = 0;	// ADD-BY-LEETEN 12/25/2012

			/////////////////////////////////////////////////////////////////
			vector<size_t> vuPoolDimLengths;
			vuPoolDimLengths.resize(this->UGetNrOfDims());
			this->vpcCoefPools.resize(this->uNrOfUpdatingCoefs);	// allocate the pools
			for(size_t c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				vector<size_t> vuPoolSubs;
				_ConvertIndexToSub(c, vuPoolSubs, this->vuDimLevels);

				bool bIsOutOfCore = false;

				// decide the pool size
				for(size_t d = 0; d < this->UGetNrOfDims(); d++)
				{
					vuPoolDimLengths[d] = (!vuPoolSubs[d])?1:(1<<vuPoolSubs[d] - 1);
					// decide whether the array is sparse.
					// MOD-BY-LEETEN 12/25/2012-FROM:	if( vuPoolSubs[d] > vuOptimalDimLevel[d] )
					if( vuPoolSubs[d] >= vuOptimalDimLevel[d] )
					// MOD-BY-LEETEN 12/25/2012-END
						bIsOutOfCore  = true;
				}

				// ADD-BY-LEETEN 12/25/2012-BEGIN
				if(bIsOutOfCore )
					uNrOfOutOfCores++;
				// ADD-BY-LEETEN 12/25/2012-END

				if( !uNrOfCoefsInFullArray )
					bIsOutOfCore  = true;

				if( !bIsOutOfCore )
				{
					vpcCoefPools[c] = new CSepDWTPool<T, unsigned short>;
					this->vpcCoefPools[c]->_Set(
						UGetNrOfBins(),
						vuPoolDimLengths,
						vuMaxCounts[c],
						true);
				} 
				else
					vpcCoefPools[c] = NULL;
			}
			LOG_VAR(uNrOfOutOfCores);	// ADD-BY-LEETEN 12/25/2012
		}

		
		// ADD-BY-LEETEN 12/23/2012-BEGIN
		//! Compute and display statistics for the computed wavelet coefficients.
		virtual	
		void
		_ShowStatistics
		(
			void *_Reserved = NULL
		)
		{
			// read the file size
			FILE* fp = fopen(szFilepath, "rb");
			fseek(fp, 0, SEEK_END);
			size_t uFileSize = ftell(fp);
			fclose(fp);
			LOG_VAR(uFileSize);

			for(size_t d = 0; d < UGetNrOfDims(); d++)
				LOG_VAR(vuDimLengths[d]);
			LOG_VAR(UGetNrOfBins());

			double dOverhead = (double)uFileSize / (double)this->uDataSize;
			LOG_VAR(dOverhead);

			double dCR = (double)(this->uDataSize * UGetNrOfBins() * sizeof(T)) / (double)uFileSize;
			LOG_VAR(dCR);

			_ShowMemoryUsage(false);
		}
		// ADD-BY-LEETEN 12/23/2012-END

		virtual 
		void
		_LoadFile
		(
			const char* szFilepath,
			void *_Reserved = NULL
		)
		{
			this->szFilepath = szFilepath;			// ADD-BY-LEETEN 12/23/2012

			/////////////////////////////////////////////////////////////////
			// now load the coefficients

			// #if WITH_NETCDF 

			#if !WITH_NETCDF4 
			ASSERT_NETCDF(nc_open(
    				szFilepath,
    				NC_NOWRITE,
    				&iNcId));
			#else	// #if !WITH_NETCDF4
			ASSERT_NETCDF(nc_open(
    				szFilepath,
    				NC_NOWRITE | NC_NETCDF4,
    				&iNcId));
			#endif	// #if !WITH_NETCDF4

			size_t uNrOfDims;
			ASSERT_NETCDF(nc_inq_dimid (
				iNcId,
				szDimDim,
				&ncDimDim));
			ASSERT_NETCDF(nc_inq_dimlen(
				iNcId,
				ncDimDim,
				&uNrOfDims));

			size_t uNrOfBins;
			ASSERT_NETCDF(nc_inq_dimid (
				iNcId,
				szDimBin,
				&ncDimBin));
			ASSERT_NETCDF(nc_inq_dimlen(
				iNcId,
				ncDimBin,
				&uNrOfBins));

			ASSERT_NETCDF(nc_inq_dimid (
				iNcId,
				szDimValue,
				&ncDimValue));

			vector<size_t> vuDimLengths;
			for(size_t t = 0; t < NR_OF_DIM_TYPES; t++)
				for(size_t d = 0; d < uNrOfDims; d++)
				{
					char szDimName[NC_MAX_NAME+1];
					sprintf(szDimName, "%s_DIM_%d", pszDimTypes[t], (unsigned int)d);
					int ncDim;
					ASSERT_NETCDF(nc_inq_dimid (
						iNcId,
						szDimName,
						&ncDim));
					vncDims.push_back(ncDim);

					size_t uDimLength;
					ASSERT_NETCDF(nc_inq_dimlen(
						iNcId,
						ncDim,
						&uDimLength));

					switch(t)
					{
					case DIM_TYPE_DATA:
						vuDimLengths.push_back(uDimLength);
						break;
					case DIM_TYPE_COEF:
						piCoefDimIds[uNrOfDims - 1 - d] = ncDim;
						break;
					}
				}

			_Set(vuDimLengths, uNrOfBins);
			_Allocate();

			//
			// now define the variable for the coef headers
			ASSERT_NETCDF(nc_inq_varid(
				iNcId,
				szVarHeaderCount,
				&ncVarHeaderCount));
			ASSERT_NETCDF(nc_inq_varid(
				iNcId,
				szVarHeaderOffset,
				&ncVarHeaderOffset));
			ASSERT_NETCDF(nc_inq_varid(
				iNcId,
				szVarCoefValue,
				&ncVarCoefValue));
			ASSERT_NETCDF(nc_inq_varid(
				iNcId,
				szVarCoefBin,
				&ncVarCoefBin));

			/////////////////////////////////////////////////////////////////////
			// define the #non-zero bins per coef.
			TBuffer<TYPE_HEADER_OFFSET> pHeaderOffset;
			TBuffer<TYPE_HEADER_COUNT> pHeaderCount;

			// now read the entire header
			pHeaderOffset.alloc(vuHeaderOffset.size());
			ASSERT_NETCDF(nc_get_var(
				iNcId,
				ncVarHeaderOffset,
				(void*)&pHeaderOffset[0]));
			for(size_t h = 0; h < vuHeaderOffset.size(); h++)
				vuHeaderOffset[h] = (size_t)pHeaderOffset[h];
			pHeaderOffset.free();

			pHeaderCount.alloc(vusHeaderCount.size());
			ASSERT_NETCDF(nc_get_var(
				iNcId,
				ncVarHeaderCount,
				(void*)&pHeaderCount[0]));
			for(size_t h = 0; h < vusHeaderCount.size(); h++)
				vusHeaderCount[h] = (unsigned short)pHeaderCount[h];
			pHeaderCount.free();

			/////////////////////////////////////////////////////////////////
			// now load the coefficients that can be in core
			vector<size_t> vuBasisHeaderSize;
			vuBasisHeaderSize.resize(UGetNrOfDims());
			vector<size_t> vuBasisHeaderBase;
			vuBasisHeaderBase.resize(UGetNrOfDims());
			size_t puStart[NC_MAX_DIMS];
			size_t puCount[NC_MAX_DIMS];
			for(size_t c = 0; c < this->uNrOfUpdatingCoefs; c++)
			{
				if( !this->vpcCoefPools[c] )
					continue;

				vector<size_t> vuLevelSub;
				_ConvertIndexToSub(c, vuLevelSub, this->vuDimMaxLevels);

				size_t uNrOfBasisHeader = 1;
				for(size_t d = 0; d < vuLevelSub.size(); d++)
				{
					// From this subscript, we can get the dim. length in this basis. 
					size_t uLevel = vuLevelSub[d];
					size_t uLen = (!uLevel)?1:(1<<(uLevel - 1));

					vuBasisHeaderSize[d] = uLen;
					uNrOfBasisHeader *= uLen;

					// we can also decide it base in the n-dim. pool
					vuBasisHeaderBase[d] = (!uLevel)?0:uLen;

					puStart[UGetNrOfDims() - 1 - d] = vuBasisHeaderBase[d];
					puCount[UGetNrOfDims() - 1 - d] = vuBasisHeaderSize[d];
				}
				pHeaderOffset.alloc(uNrOfBasisHeader);
				pHeaderCount.alloc(uNrOfBasisHeader);

				// read the header for the current basis
				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarHeaderOffset,
					puStart,
					puCount,
					(void*)&pHeaderOffset[0]));

				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarHeaderCount,
					puStart,
					puCount,
					(void*)&pHeaderCount[0]));

				// compute the total number of coefficients in the current basis
				size_t uNrOfBasisCoefs = 0;
				for(size_t h = 0; h < uNrOfBasisHeader; h++)
					uNrOfBasisCoefs += (size_t)pHeaderCount[h];
				puStart[0] = (size_t)pHeaderOffset[0];
				puCount[0] = uNrOfBasisCoefs;

				// read the coefficients of the current basis
				TBuffer<TYPE_COEF_BIN> pBasisCoefBin;
				pBasisCoefBin.alloc(uNrOfBasisCoefs);
				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarCoefBin,
					puStart,
					puCount,
					(void*)&pBasisCoefBin[0]));

				TBuffer<TYPE_COEF_VALUE> pBasisCoefValue;
				pBasisCoefValue.alloc(uNrOfBasisCoefs);
				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarCoefValue,
					puStart,
					puCount,
					(void*)&pBasisCoefValue[0]));

				// scan through all basis
				for(size_t coefi = 0, basis = 0; basis < uNrOfBasisHeader; basis++)
				{
					vector<size_t> vuBasisSub;
					_ConvertIndexToSub(basis, vuBasisSub, vuBasisHeaderSize);

					// scan through all bin
					for(size_t bini = 0; bini < pHeaderCount[basis]; bini++, coefi++)
						this->vpcCoefPools[c]->_AddAt(
							pBasisCoefBin[coefi],
							vuBasisSub,
							pBasisCoefValue[coefi]);
				}
				this->vpcCoefPools[c]->_Finalize(1.0);
			}
			// #endif	// #if WITH_NETCDF 
		}

		//! Return the sum of all bins at the given position
		virtual	
		void
		_GetAllSums
		(
			const vector<size_t>& vuPos,
			vector<double>& vdSums,
			void *_Reserved = NULL
		)
		{
			vector<size_t> vuEmpty;	// this empty vector is used to speed up the access of elements

			vdSums.resize(UGetNrOfBins());

			#if	0	// DEL-BY-LEETEN 12/25/2012-BEGIN
			vector<double> vdWaveletBasis;
			vdWaveletBasis.resize( this->uNrOfWaveletsFromAllDims );

			vector<size_t> vuSubs;		vuSubs.resize( this->uNrOfWaveletsFromAllDims );
			#endif	// DEL-BY-LEETEN 12/25/2012-END
			vector<size_t> vuPoolSub;	vuPoolSub.resize( UGetNrOfDims() );
			vector<size_t> vuPoolBase;	vuPoolBase.resize( UGetNrOfDims() );
			vector<size_t> vuSubInPool;	vuSubInPool.resize( UGetNrOfDims() );
			vector<size_t> vuHeaderSub;	vuHeaderSub.resize( UGetNrOfDims() );
		
			#if	0	// MOD-BY-LEETEN 12/25/2012-FROM:
			// for each dimenion d, based on the posistion, store the corresponding l[d] wavelet basis value
			for(size_t p = 0, d = 0; d < UGetNrOfDims(); d++)
			{
				size_t uPos = vuPos[d];
				size_t uDimMaxLevel = vuDimMaxLevels[d];
				size_t uMaxWin = 1 << vuDimLevels[d];
				for(size_t 	
					l = 0, w = uMaxWin; 
					l < uDimMaxLevel; 
					l++, p++, w >>= 1)
				{
					// Decide the wavelet size based on the current level, and then the wavelet basis values based on the position within this wavelet
					double dWaveletBasis = 0.0;
					
					// assume Haar wavelet for now
					// Given the currenet level l and subscript uPos, compute the corresponding wavelet value
					size_t uPosInWavelet = uPos % w;				
					if( uPosInWavelet < w / 2)
						dWaveletBasis = 1.0;
					else
						dWaveletBasis = -1.0;
						
					if( l >= 2 )
						dWaveletBasis *= sqrt( (double)(1 << (l - 1)) );
					vdWaveletBasis[p] = dWaveletBasis;

					vuSubs[p] = uPos / w;
				}
			}
			#else	// MOD-BY-LEETEN 12/25/2012-TO:
			vector<double> vdWaveletBasis;
			vector<size_t> vuSubs;
			_GetBackwardWavelet(
				vuPos, 
				vuSubs, 
				vdWaveletBasis, 
				true);
			#endif	// MOD-BY-LEETEN 12/25/2012-END

			// now find the combination of the coefficients of all dimensions 
			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				_ConvertIndexToSub(c, vuPoolSub, this->vuDimLevels);

				double dWavelet = 1.0;
				for(size_t d = 0, uBase = 0;
					d < UGetNrOfDims(); 
					uBase += vuDimMaxLevels[d], d++, p++)
				{
					// compute the wavelet
					dWavelet *= vdWaveletBasis[uBase + vuCoefDim2Level[p]];	

					// compute the subscript and size for the current basis
					vuSubInPool[d] = vuSubs[uBase + vuCoefDim2Level[p]];

					// 
					vuPoolBase[d] = (!vuPoolSub[d])?0:(1<<(vuPoolSub[d]-1));

					// decide the subscript in the 
					vuHeaderSub[d] = vuPoolBase[d] + vuSubInPool[d];
				}

				if( this->vpcCoefPools[c] )
				{
					vector< pair<size_t, double> > vpairCoefs;
					this->vpcCoefPools[c]->_GetCoefSparse
					(
						vuSubInPool,
						vpairCoefs
					);

					for(vector< pair<size_t, double> >::iterator
						ivpairCoefs = vpairCoefs.begin();
						ivpairCoefs != vpairCoefs.end();
						ivpairCoefs++ )
						vdSums[ivpairCoefs->first] += ivpairCoefs->second * dWavelet;
				}
				else
				{
					size_t uHeaderIndex = UConvertSubToIndex(vuHeaderSub, this->vuCoefLengths);
					size_t uStart = (size_t)vuHeaderOffset[uHeaderIndex];
					size_t uCount = (size_t)vusHeaderCount[uHeaderIndex];

					// ADD-BY-LEETEN 12/25/2012-BEGIN
					if( uCount )
					{
					// ADD-BY-LEETEN 12/25/2012-END

					ASSERT_NETCDF(nc_get_vara(
						iNcId,
						ncVarCoefBin,
						&uStart,
						&uCount,
						(void*)&pCoefBin[0]));
	
					ASSERT_NETCDF(nc_get_vara(
						iNcId,
						ncVarCoefValue,
						&uStart,
						&uCount,
						(void*)&pCoefValue[0]));
					}	// ADD-BY-LEETEN 12/25/2012

					for(size_t i = 0; i < uCount; i++)
						vdSums[pCoefBin[i]] += pCoefValue[i] * dWavelet;
				}
			}
			for(size_t b = 0; b < UGetNrOfBins(); b++)
				vdSums[b] /= dWaveletDenomiator;
		}

		CSATSepDWTOutOfCore()
		{
		}

		~CSATSepDWTOutOfCore()
		{
			for(size_t c = 0; c < uNrOfUpdatingCoefs; c++)
					if(vpcCoefPools[c])
						delete vpcCoefPools[c];

			/*
			// Ideally, the file should be closed. Nevertheless, I will get an error message at this point. 
			As this file is read only, it should be fine NOT to close it.
			ASSERT_NETCDF(nc_close(iNcId));
			*/
		}
	};
}
