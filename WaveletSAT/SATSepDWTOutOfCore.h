#pragma once

// ADD-BY-LEETEN 12/28/2012-BEGIN
#if	WITH_BOOST
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;
#endif	// #if	WITH_BOOST
// ADD-BY-LEETEN 12/28/2012-END

#include <map>	

#include <vector>
#include <algorithm>	// ADD-BY-LEETEN 12/26/2012
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

// ADD-BY-LEETEN 12/28/2012-BEGIN
// if this is non-0, the header with most coefficients will be place in core first
// MOD-BY-LEETEN 12/29/2012-FROM:	#define	IS_SELECTING_MOST_FIRST						1
#define	IS_SELECTING_LONGEST_FIRST					1
// MOD-BY-LEETEN 12/29/2012-END

// if this is non-0, when a coefficient is chosen, other coefficients of the same wavelet will be chosen as well
#define IS_SELECTING_THE_SAME_WAVELET			1
// ADD-BY-LEETEN 12/28/2012-END

namespace WaveletSAT
{
	// MOD-BY-LEETEN 12/29/2012-FROM:	template<typename T>
	template<
		typename DT,	//!< Type of the data
		typename WT		//!< Type of the wavelet coefficients
	>
	// MOD-BY-LEETEN 12/29/2012-END
	// The class that access the coefficients from files (in NetCDF format)
	class CSATSepDWTOutOfCore:
		public CSATSepDWTNetCDF,
		public CSepDWTHeader
	{
protected:	
			// ADD-BY-LEETEN 12/28/2012-BEGIN
			//! The accumulated #I/O requests since the I/O counters are reset.
			size_t uAccumNrOfIORequest;

			//! The max #I/O requests per query since the I/O counters are reset.
			size_t uMaxNrOfIORequest;

			//! The min #I/O requests per query since the I/O counters are reset.
			size_t uMinNrOfIORequest;

			//! The #query since the I/O counters are reset.
			size_t uNrOfQueries;
			// ADD-BY-LEETEN 12/28/2012-END

			const char* szFilepath;	

			#if	0	// MOD-BY-LEETEN 12/29/2012-FROM:
			vector<size_t>			vuHeaderOffset;
			vector<unsigned short>	vusHeaderCount;
			#else	// MOD-BY-LEETEN 12/29/2012-TO:
			//! The D-dim. array of the offset to the 1D array of coefficients
			vector<size_t>			vuCoefOffsets;

			//! The D-dim. array of the coefficient counters
			vector<unsigned short>	vusCoefCounts;
			#endif	// MOD-BY-LEETEN 12/29/2012-END

			#if	0	// MOD-BY-LEETEN 12/29/2012-FROM:
			//! #Coefs stored in full arrays
			size_t uNrOfCoefsInFullArray;	
			#else	// MOD-BY-LEETEN 12/29/2012-TO:
			//! #Coef. values stored in full arrays
			size_t uMaxNrOfValuesInCore;	
			#endif	// MOD-BY-LEETEN 12/29/2012-END

			//! vector of pointers to the coefficient pools. 
			/*!
			If the poiner is NULL, it means that the coefficients are out of core.
			*/
			vector< CSepDWTPool<WT, unsigned short>* > vpcCoefPools;

			#if	0	// MOD-BY-LEETEN 12/29/2012-FROM:
			TBuffer<TYPE_COEF_BIN>		pCoefBin;
			TBuffer<TYPE_COEF_VALUE>	pCoefValue;

			// ADD-BY-LEETEN 12/26/2012-BEGIN
			size_t uNrOfNonZeroCoefs;
			vector<bool> vbHeaderInCore;
			#else	// MOD-BY-LEETEN 12/29/2012-TO:
			//! An array to store the coefficient bins.
			TBuffer<TYPE_COEF_BIN>		pCoefBins;

			//! An array to store the coefficient values.
			TBuffer<TYPE_COEF_VALUE>	pCoefValues;

			//! #Non-zero values in the file
			size_t uNrOfNonZeroValues;

			//! The D-dim array of flags indicating whether the corresponding coefficients are in core
			vector<bool> vbFlagsCoefInCore;
			#endif	// MOD-BY-LEETEN 12/29/2012-END

			// ADD-BY-LEETEN 12/26/2012-END

public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0a00,
			// ADD-BY-LEETEN 12/28/2012-BEGIN
			ACCUM_NR_OF_IO_REQUESTS,
			MAX_NR_OF_IO_REQUESTS,
			MIN_NR_OF_IO_REQUESTS,
			RESET_IO_COUNTERS,
			// ADD-BY-LEETEN 12/28/2012-END
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

		// ADD-BY-LEETEN 12/28/2012-BEGIN
		virtual	
		void
		_GetLong(
			int eName,
			long *plValue,
			void* _Reserved = NULL
		)
		{
			switch(eName)
			{
			case ACCUM_NR_OF_IO_REQUESTS:
				*plValue = (long)uAccumNrOfIORequest;
				break;
			case MAX_NR_OF_IO_REQUESTS:
				*plValue = (long)uMaxNrOfIORequest;
				break;
			case MIN_NR_OF_IO_REQUESTS:
				*plValue = (long)uMinNrOfIORequest;
				break;
			}
		}

		virtual
		void
		_SetBoolean
		(
			int eName,
			bool bValue,
			void* _Reserved = NULL
		)
		{
			switch(eName)
			{
			case RESET_IO_COUNTERS:
				uAccumNrOfIORequest = 0;
				uMaxNrOfIORequest = 0;
				uMinNrOfIORequest = uNrOfUpdatingCoefs;
				uNrOfQueries = 0;
				break;
			}
		}
		// ADD-BY-LEETEN 12/28/2012-END

		virtual
		void
		_Allocate(
					void *_Reserved = NULL
					)
		{
			#if	0	// DEL-BY-LEETEN 12/26/2012-BEGIN
			vuCoefOffsets.resize(uNrOfCoefs);
			vusCoefCounts.resize(uNrOfCoefs);
			#endif	// DEL-BY-LEETEN 12/26/2012-END
			pCoefBins.alloc(UGetNrOfBins());
			pCoefValues.alloc(UGetNrOfBins());

			#if	0	// MOD-BY-LEETEN 12/26/2012-FROM:
			// decide which basis should be in core
			uNrOfValuesInCore = min(
				(size_t)floor( (double)uSizeOfFullArrays/(double)(UGetNrOfBins() * sizeof(WT))), 
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
			
			if( uNrOfValuesInCore )
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

					if( uSize <= uNrOfValuesInCore )
					{
						size_t uNewDiff = uNrOfValuesInCore - uSize;
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

				if( !uNrOfValuesInCore )
					bIsOutOfCore  = true;

				if( !bIsOutOfCore )
				{
					vpcCoefPools[c] = new CSepDWTPool<WT, unsigned short>;
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
			#else	// MOD-BY-LEETEN 12/26/2012-TO:
			// sort the coef. indices by its # coefficients
			// MOD-BY-LEETEN 12/29/2012-FROM:	vector< pair<unsigned short, long long> > vpairHeader;
			vector< pair<unsigned short, long long> > vpairCoefCountIndex;
			// MOD-BY-LEETEN 12/29/2012-END
			vpairCoefCountIndex.resize(uNrOfCoefs);
			for(size_t c = 0; c < uNrOfCoefs; c++)
				vpairCoefCountIndex[c] = make_pair(this->vusCoefCounts[c], -(long long)c);
			sort(vpairCoefCountIndex.begin(), vpairCoefCountIndex.end());

			////////////////////////////////////////////////////////
			// now collect the coefficients
			// Mark the selected header
			uMaxNrOfValuesInCore = min(
				(size_t)floor( (double)uSizeOfFullArrays/(double)sizeof(WT)), 
				uNrOfNonZeroValues);
			vbFlagsCoefInCore.resize(uNrOfCoefs);
			#if	0	// MOD-BY-LEETEN 12/28/2012-FROM:
			for(size_t c = 0, uNrOfValuesInCore = 0; c < uNrOfCoefs && uNrOfValuesInCore < uNrOfValuesInCore; uNrOfValuesInCore += vpairCoefCountIndex[c].first, c++)
				vbFlagsCoefInCore[abs(vpairCoefCountIndex[c].second)] = true;
			#else	// MOD-BY-LEETEN 12/28/2012-TO:
			// vector<size_t> vuCoefSub;	vuCoefSub.resize(UGetNrOfDims());
			vector<size_t> vuLevel;		vuLevel.resize(UGetNrOfDims());
			vector<size_t> vuLevelBase;	vuLevelBase.resize(UGetNrOfDims());
			vector<size_t> vuLevelSize;	vuLevelSize.resize(UGetNrOfDims());
			vector<size_t> vuSubInLevel;vuSubInLevel.resize(UGetNrOfDims());
			vector<size_t> vuCoefSub;	vuCoefSub.resize(UGetNrOfDims());
			// MOD-BY-LEETEN 12/29/2012-FROM:	size_t uCollectedCoefs = 0;
			size_t uNrOfValuesInCore = 0;
			// MOD-BY-LEETEN 12/29/2012-END

			#if				IS_SELECTING_LONGEST_FIRST
			reverse(vpairCoefCountIndex.begin(), vpairCoefCountIndex.end());
			#endif	// #if	IS_SELECTING_LONGEST_FIRST

			for(size_t c = 0; c < uNrOfCoefs && uNrOfValuesInCore <= uMaxNrOfValuesInCore; c++)
			{
				#if		!IS_SELECTING_THE_SAME_WAVELET
				if( uNrOfValuesInCore + vpairCoefCountIndex[c].first > uNrOfValuesInCore )
					break;

				vbFlagsCoefInCore[abs(vpairCoefCountIndex[c].second)] = true;
				uNrOfValuesInCore += vpairCoefCountIndex[c].first; 
				#else	// #if		!IS_SELECTING_THE_SAME_WAVELET
				// MOD-BY-LEETEN 12/29/2012-FROM:	size_t uHeader = abs(vpairCoefCountIndex[c].second);
				size_t uCoef = abs(vpairCoefCountIndex[c].second);
				// MOD-BY-LEETEN 12/29/2012-END
				if( vbFlagsCoefInCore[uCoef] )
					continue;

				_ConvertIndexToLevels
				(
					uCoef,
					vuLevel,
					vuSubInLevel,
					vuLevelBase,
					vuLevelSize
				);
				size_t uLevelSize = 1;
				for(size_t d = 0; d < vuLevelSize.size(); d++)
					uLevelSize *= vuLevelSize[d];

				for(size_t i = 0; i < uLevelSize && uNrOfValuesInCore < uMaxNrOfValuesInCore; i++)
				{
					_ConvertIndexToSub(i, vuSubInLevel, vuLevelSize);
					for(size_t d = 0; d < UGetNrOfDims(); d++)
						vuCoefSub[d] = vuLevelBase[d] + vuSubInLevel[d];
					size_t uCoefIndex = UConvertSubToIndex(vuCoefSub, vuCoefLengths);
					if( uNrOfValuesInCore + vusCoefCounts[uCoefIndex] > uMaxNrOfValuesInCore )
						break;
					vbFlagsCoefInCore[uCoefIndex] = true;
					uNrOfValuesInCore += vusCoefCounts[uCoefIndex];
				}
				#endif	// #if		!IS_SELECTING_THE_SAME_WAVELET
			}
			LOG_VAR(uNrOfNonZeroValues);
			LOG_VAR(uMaxNrOfValuesInCore);
			LOG_VAR(uNrOfValuesInCore);
			#endif	// MOD-BY-LEETEN 12/28/2012-END

			this->vpcCoefPools.resize(this->uNrOfUpdatingCoefs);	// allocate the pools
			#endif	// MOD-BY-LEETEN 12/26/2012-END
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
			// ADD-BY-LEETEN 12/28/2012-BEGIN
			#if		WITH_BOOST
			fs::path pathNativePath( szFilepath, fs::native );
			size_t uFileSize = fs::file_size( pathNativePath );
			#else	// #if WITH_BOOST
			// ADD-BY-LEETEN 12/28/2012-END
			FILE* fp = fopen(szFilepath, "rb");
			fseek(fp, 0, SEEK_END);
			size_t uFileSize = ftell(fp);
			fclose(fp);
			#endif	// #if WITH_BOOST	// ADD-BY-LEETEN 12/28/2012
			LOG_VAR(uFileSize);

			for(size_t d = 0; d < UGetNrOfDims(); d++)
				LOG_VAR(vuDimLengths[d]);
			LOG_VAR(UGetNrOfBins());

			double dOverhead = (double)uFileSize / (double)this->uDataSize;
			LOG_VAR(dOverhead);

			double dCR = (double)(this->uDataSize * UGetNrOfBins() * sizeof(WT)) / (double)uFileSize;
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
			// ADD-BY-LEETEN 12/26/2012-BEGIN
			ASSERT_NETCDF(nc_inq_dimlen (
				iNcId,
				ncDimValue,
				&uNrOfNonZeroValues));
			// ADD-BY-LEETEN 12/26/2012-END

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
			// DEL-BY-LEETEN 12/26/2012:	_Allocate();

			//
			// now define the variable for the coef headers
			ASSERT_NETCDF(nc_inq_varid(
				iNcId,
				szVarCoefCount,
				&ncVarCoefCount));
			ASSERT_NETCDF(nc_inq_varid(
				iNcId,
				szVarCoefOffset,
				&ncVarCoefOffset));
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
			#if	0	// MOD-BY-LEETEN 12/29/2012-FROM:
			TBuffer<TYPE_COEF_OFFSET> pHeaderOffset;
			TBuffer<TYPE_COEF_COUNT> pHeaderCount;
			#else	// MOD-BY-LEETEN 12/29/2012-TO:
			TBuffer<TYPE_COEF_OFFSET> pCoefOffsets;
			TBuffer<TYPE_COEF_COUNT> pCoefCounts;
			#endif	// MOD-BY-LEETEN 12/29/2012-END

			// now read the entire header
			// ADD-BY-LEETEN 12/26/2012-BEGIN
			vuCoefOffsets.resize(uNrOfCoefs);
			vusCoefCounts.resize(uNrOfCoefs);
			// ADD-BY-LEETEN 12/26/2012-END
			pCoefOffsets.alloc(vuCoefOffsets.size());
			ASSERT_NETCDF(nc_get_var(
				iNcId,
				ncVarCoefOffset,
				(void*)&pCoefOffsets[0]));
			for(size_t h = 0; h < vuCoefOffsets.size(); h++)
				vuCoefOffsets[h] = (size_t)pCoefOffsets[h];
			pCoefOffsets.free();

			pCoefCounts.alloc(vusCoefCounts.size());
			ASSERT_NETCDF(nc_get_var(
				iNcId,
				ncVarCoefCount,
				(void*)&pCoefCounts[0]));
			for(size_t h = 0; h < vusCoefCounts.size(); h++)
				vusCoefCounts[h] = (unsigned short)pCoefCounts[h];
			pCoefCounts.free();

			// ADD-BY-LEETEN 12/26/2012-BEGIN
			_Allocate();
			size_t uNrOfAllocatedPools = 0;
			// DEL-BY-LEETEN 12/29/2012:	size_t uNrOfInCoreHeaders = 0;
			// ADD-BY-LEETEN 12/26/2012-END

			/////////////////////////////////////////////////////////////////
			// now load the coefficients that can be in core
			// MOD-BY-LEETEN 12/29/2012-FROM:	vector<size_t> vuBasisHeaderSize;
			vector<size_t> vuLocalCoefLengths;
			// MOD-BY-LEETEN 12/29/2012-END
			vuLocalCoefLengths.resize(UGetNrOfDims());
			// MOD-BY-LEETEN 12/29/2012-FROM:	vector<size_t> vuBasisHeaderBase;
			vector<size_t> vuGlobalCoefBase;
			// MOD-BY-LEETEN 12/29/2012-END
			vuGlobalCoefBase.resize(UGetNrOfDims());

			size_t puStart[NC_MAX_DIMS];
			size_t puCount[NC_MAX_DIMS];
			for(size_t c = 0; c < this->uNrOfUpdatingCoefs; c++)
			{
				#if	0	// DEL-BY-LEETEN 12/26/2012-BEGIN
				if( !this->vpcCoefPools[c] )
					continue;
				#endif	// DEL-BY-LEETEN 12/26/2012-END

				vector<size_t> vuLevelSub;
				_ConvertIndexToSub(c, vuLevelSub, this->vuDimMaxLevels);

				// MOD-BY-LEETEN 12/29/2012-FROM:				size_t uNrOfBasisHeader = 1;
				size_t uNrOfLocalCoefs = 1;
				// MOD-BY-LEETEN 12/29/2012-END
				for(size_t d = 0; d < vuLevelSub.size(); d++)
				{
					// From this subscript, we can get the dim. length in this basis. 
					size_t uLevel = vuLevelSub[d];
					size_t uLen = (!uLevel)?1:(1<<(uLevel - 1));

					vuLocalCoefLengths[d] = uLen;
					uNrOfLocalCoefs *= uLen;

					// we can also decide it base in the n-dim. pool
					vuGlobalCoefBase[d] = (!uLevel)?0:uLen;

					puStart[UGetNrOfDims() - 1 - d] = vuGlobalCoefBase[d];
					puCount[UGetNrOfDims() - 1 - d] = vuLocalCoefLengths[d];
				}
				// ADD-BY-LEETEN 12/29/2012-BEGIN
				TBuffer<TYPE_COEF_OFFSET> pLocalCoefOffsets;
				TBuffer<TYPE_COEF_COUNT> pLocalCoefCounts;
				// ADD-BY-LEETEN 12/29/2012-END
				pLocalCoefOffsets.alloc(uNrOfLocalCoefs);
				pLocalCoefCounts.alloc(uNrOfLocalCoefs);

				// read the header for the current basis
				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarCoefOffset,
					puStart,
					puCount,
					(void*)&pLocalCoefOffsets[0]));

				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarCoefCount,
					puStart,
					puCount,
					(void*)&pLocalCoefCounts[0]));

				// compute the total number of coefficients in the current basis
				// MOD-BY-LEETEN 12/29/2012-FROM:	size_t uNrOfBasisCoefs = 0;
				size_t uNrOfLocalCoefValues = 0;
				// MOD-BY-LEETEN 12/29/2012-END
				for(size_t h = 0; h < uNrOfLocalCoefs; h++)
					uNrOfLocalCoefValues += (size_t)pLocalCoefCounts[h];
				puStart[0] = (size_t)pLocalCoefOffsets[0];
				puCount[0] = uNrOfLocalCoefValues;

				// read the coefficients of the current basis
				// MOD-BY-LEETEN 12/29/2012-FROM:	TBuffer<TYPE_COEF_BIN> pBasisCoefBin;
				TBuffer<TYPE_COEF_BIN> pLocalCoefBins;
				// MOD-BY-LEETEN 12/29/2012-END
				pLocalCoefBins.alloc(uNrOfLocalCoefValues);
				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarCoefBin,
					puStart,
					puCount,
					(void*)&pLocalCoefBins[0]));

				// MOD-BY-LEETEN 12/29/2012-FROM:				TBuffer<TYPE_COEF_VALUE> pBasisCoefValue;
				TBuffer<TYPE_COEF_VALUE> pLocalCoefValues;
				// MOD-BY-LEETEN 12/29/2012-END
				pLocalCoefValues.alloc(uNrOfLocalCoefValues);
				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarCoefValue,
					puStart,
					puCount,
					(void*)&pLocalCoefValues[0]));

				// scan through all basis
				// MOD-BY-LEETEN 12/29/2012-FROM:				for(size_t coefi = 0, basis = 0; basis < uNrOfLocalCoefs; basis++)
				for(size_t valuei = 0, lc = 0; lc < uNrOfLocalCoefs; lc++)
				// MOD-BY-LEETEN 12/29/2012-END
				{
					// MOD-BY-LEETEN 12/29/2012-FROM:			vector<size_t> vuBasisSub;
					vector<size_t> vuLocalCoefSub;
					// MOD-BY-LEETEN 12/29/2012-END
					_ConvertIndexToSub(lc, vuLocalCoefSub, vuLocalCoefLengths);

					// ADD-BY-LEETEN 12/26/2012-BEGIN
					vector<size_t> vuCoefSub;
					vuCoefSub.resize(UGetNrOfDims());
					for(size_t d = 0; d < UGetNrOfDims(); d++)
						vuCoefSub[d] = vuGlobalCoefBase[d] + vuLocalCoefSub[d];
					size_t uCoefIndex = UConvertSubToIndex(vuCoefSub, this->vuCoefLengths);
					bool bIsInCore = this->vbFlagsCoefInCore[uCoefIndex];

					if( !bIsInCore )
					{
						valuei += pLocalCoefCounts[lc];
						continue;
					}

					// DEL-BY-LEETEN 12/29/2012:	uNrOfInCoreHeaders++;
					if( !vpcCoefPools[c] )
					{
						vpcCoefPools[c]= new CSepDWTPool<WT, unsigned short>;
						vpcCoefPools[c]->_Set(
							UGetNrOfBins(),
							vuLocalCoefLengths,
							vuMaxCounts[c],
							true);
						uNrOfAllocatedPools++;
					}
					// ADD-BY-LEETEN 12/26/2012-END

					// scan through all bin
					for(size_t bini = 0; bini < pLocalCoefCounts[lc]; bini++, valuei++)
						this->vpcCoefPools[c]->_AddAt(
							pLocalCoefBins[valuei],
							vuLocalCoefSub,
							pLocalCoefValues[valuei]);
				}
				if( vpcCoefPools[c] )	// ADD-BY-LEETEN 12/26/2012
				this->vpcCoefPools[c]->_Finalize(1.0);
			}
			// ADD-BY-LEETEN 12/26/2012-BEGIN
			#if	0	// DEL-BY-LEETEN 12/29/2012-BEGIN
			LOG_VAR(uNrOfCoefs);
			LOG_VAR(uNrOfInCoreHeaders);
			#endif	// DEL-BY-LEETEN 12/29/2012-END
			LOG_VAR(uNrOfAllocatedPools);
			// ADD-BY-LEETEN 12/26/2012-END
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
			#if	0	// MOD-BY-LEETEN 12/29/2012-FROM:
			vector<size_t> vuPoolSub;	vuPoolSub.resize( UGetNrOfDims() );
			vector<size_t> vuPoolBase;	vuPoolBase.resize( UGetNrOfDims() );
			vector<size_t> vuSubInPool;	vuSubInPool.resize( UGetNrOfDims() );
			vector<size_t> vuHeaderSub;	vuHeaderSub.resize( UGetNrOfDims() );
			#else	// MOD-BY-LEETEN 12/29/2012-TO:
			vector<size_t> vuLevels;			vuLevels.resize( UGetNrOfDims() );
			vector<size_t> vuGlobalCoefBase;	vuGlobalCoefBase.resize( UGetNrOfDims() );
			vector<size_t> vuLocalCoefSub;		vuLocalCoefSub.resize( UGetNrOfDims() );
			vector<size_t> vuCoefSub;			vuCoefSub.resize( UGetNrOfDims() );
			#endif	// MOD-BY-LEETEN 12/29/2012-END
		
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

			// ADD-BY-LEETEN 12/28/2012-BEGIN
			uNrOfQueries++;
			size_t uNrOfIORequest = 0;
			// ADD-BY-LEETEN 12/28/2012-END

			// now find the combination of the coefficients of all dimensions 
			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				_ConvertIndexToSub(c, vuLevels, this->vuDimLevels);

				double dWavelet = 1.0;
				for(size_t d = 0, uBase = 0;
					d < UGetNrOfDims(); 
					uBase += vuDimMaxLevels[d], d++, p++)
				{
					// compute the wavelet
					dWavelet *= vdWaveletBasis[uBase + vuCoefDim2Level[p]];	

					// compute the subscript and size for the current basis
					vuLocalCoefSub[d] = vuSubs[uBase + vuCoefDim2Level[p]];

					// 
					vuGlobalCoefBase[d] = (!vuLevels[d])?0:(1<<(vuLevels[d]-1));

					// decide the subscript in the 
					vuCoefSub[d] = vuGlobalCoefBase[d] + vuLocalCoefSub[d];
				}

				// MOD-BY-LEETEN 12/26/2012-FROM:	if( this->vpcCoefPools[c] )
				// MOD-BY-LEETEN 12/29/2012-FROM:	size_t uHeaderIndex = UConvertSubToIndex(vuCoefSub, this->vuCoefLengths);
				size_t uGlobalCoefIndex = UConvertSubToIndex(vuCoefSub, this->vuCoefLengths);
				// MOD-BY-LEETEN 12/29/2012-END
				if( this->vbFlagsCoefInCore[uGlobalCoefIndex] )
				// MOD-BY-LEETEN 12/26/2012-END
				{
					// MOD-BY-LEETEN 12/29/2012-FROM:	vector< pair<size_t, double> > vpairCoefs;
					vector< pair<size_t, WT> > vpairCoefBinValues;
					// MOD-BY-LEETEN 12/29/2012-END
					this->vpcCoefPools[c]->_GetCoefSparse
					(
						vuLocalCoefSub,
						vpairCoefBinValues
					);

					for(vector< pair<size_t, double> >::iterator
						ivpairCoefs = vpairCoefBinValues.begin();
						ivpairCoefs != vpairCoefBinValues.end();
						ivpairCoefs++ )
						vdSums[ivpairCoefs->first] += ivpairCoefs->second * dWavelet;
				}
				else
				{
					uNrOfIORequest++;	// ADD-BY-LEETEN 12/28/2012

					// DEL-BY-LEETEN 12/26/2012:	size_t uGlobalCoefIndex = UConvertSubToIndex(vuCoefSub, this->vuCoefLengths);
					size_t uStart = (size_t)vuCoefOffsets[uGlobalCoefIndex];
					size_t uCount = (size_t)vusCoefCounts[uGlobalCoefIndex];

					// ADD-BY-LEETEN 12/25/2012-BEGIN
					if( uCount )
					{
					// ADD-BY-LEETEN 12/25/2012-END

					ASSERT_NETCDF(nc_get_vara(
						iNcId,
						ncVarCoefBin,
						&uStart,
						&uCount,
						(void*)&pCoefBins[0]));
	
					ASSERT_NETCDF(nc_get_vara(
						iNcId,
						ncVarCoefValue,
						&uStart,
						&uCount,
						(void*)&pCoefValues[0]));
					}	// ADD-BY-LEETEN 12/25/2012

					for(size_t i = 0; i < uCount; i++)
						vdSums[pCoefBins[i]] += pCoefValues[i] * dWavelet;
				}
			}
			// ADD-BY-LEETEN 12/28/2012-BEGIN
			uAccumNrOfIORequest += uNrOfIORequest;
			uMaxNrOfIORequest = max(uMaxNrOfIORequest, uNrOfIORequest);
			uMinNrOfIORequest = min(uMinNrOfIORequest, uNrOfIORequest);
			// ADD-BY-LEETEN 12/28/2012-END
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
					{	// ADD-BY-LEETEN 12/26/2012
						delete vpcCoefPools[c];
					// ADD-BY-LEETEN 12/26/2012-BEGIN
						vpcCoefPools[c] = NULL;
					}
					// ADD-BY-LEETEN 12/26/2012-END

			/*
			// Ideally, the file should be closed. Nevertheless, I will get an error message at this point. 
			As this file is read only, it should be fine NOT to close it.
			ASSERT_NETCDF(nc_close(iNcId));
			*/
		}
	};
}
