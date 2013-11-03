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
#include <assert.h>		// ADD-BY-LEETEN 02/03/2013

#include "SepDWTHeader.h"
#include "SATSepDWTNetCDF.h"
// MOD-BY-LEETEN 2013/10/27-FROM:	#include "SepDWTPool.h"
#include "SepDWTDecoderPool.h"
// MOD-BY-LEETEN 2013/10/27-END
#include "DecoderBase.h"	// ADD-BY-LEETEN 01/02/2013

#include "liblog.h"	

#if	WITH_NETCDF
#include <netcdf.h>
#include "lognc.h"
#endif	// #if	WITH_NETCDF

// ADD-BY-LEETEN 01/05/2013-BEGIN
// If this is non-0, the coefficient for a region query will be merged together to reduce the access of coefficients
#define MERGE_COEFFICIENTS_PER_REGION			1
// ADD-BY-LEETEN 01/05/2013-END

namespace WaveletSAT
{
	template<
		// typename DT,				//!< Type of the data
		typename ST = typeSum,		//!< Type of the sum
		typename BT = typeBin,		//!< Type of the bin
		typename WT = typeWavelet	//!< Type of the wavelet coefficientsd
	>
	// The class that access the coefficients from files (in NetCDF format)
	class CSATSepDWTDecoder:
		virtual public CDecoderBase<ST, BT>,
		virtual public CSATSepDWTNetCDF,
		virtual public CSepDWTHeader
	{
protected:	
			bool bIsPrintingTiming;	// ADD-BY-LEETEN 01/27/2013

			vector< pair<size_t, WT> > vpairLocalWaveletQueues; 

			// ADD-BY-LEETEN 01/05/2012-BEGIN
			size_t uMinNrOfCoefQueries;
			size_t uMaxNrOfCoefQueries;
			size_t uAccumNrOfCoefQueries;
			size_t uNrOfRangeQueries;
			// ADD-BY-LEETEN 01/05/2012-END

			vector<BT>	vusCachedNextOffsets;

			// vector<unsigned short>	vusCachedBins;
			// vector<WT>				vCachedValues;

			//! The base to the global pool of values per wavelet
			vector<size_t> vuGlobalValueBase;

			//! The #values per wavelet
			vector<size_t> vuLocalValueCount;
			// ADD-BY-LEETEN 12/29/2012-END

			//! The D-dim. array of the offset to the 1D array of coefficients
			vector<size_t>			vuCoefOffsets;

			//! The D-dim. array of the coefficient counters
			vector<BT>	vusCoefCounts;

			//! #Coef. values stored in full arrays
			size_t uMaxNrOfValuesInCore;	

			//! vector of pointers to the coefficient pools. 
			/*!
			If the poiner is NULL, it means that the coefficients are out of core.
			*/

			// MOD-BY-LEETEN 2013/10/27-FROM:			vector< CSepDWTPool<WT, BT>* > vpcCoefPools;
			typedef CSepDWTDecoderPool<WT, BT> CSepDWTPool;
			vector< CSepDWTPool* > vpcCoefPools;
			// MOD-BY-LEETEN 2013/10/27-END

			vector<TYPE_COEF_BIN>	vCoefBins;

			vector<TYPE_COEF_VALUE>	vCoefValues;

			//! #Non-zero values in the file
			size_t uNrOfNonZeroValues;

			//! The D-dim array of flags indicating whether the corresponding coefficients are in core
			vector<bool> vbFlagsCoefInCore;

			// ADD-BY-LEETEN 12/26/2012-END

public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0700,
			PARAMETER_END
		};

		virtual	
		void
		_SetInteger(
			int eName,
			long lValue,
			void* _Reserved = NULL
		)
		{
			CDecoderBase<ST, BT>::_SetInteger(eName, lValue);	// ADD-BY-LEETEN 01/02/2013
			CSATSepDWTNetCDF::_SetInteger(eName, lValue);
			CSepDWTHeader::_SetInteger(eName, lValue);
			// ADD-BY-LEETEN 12/30/2012-BEGIN
			switch(eName)
			{
			}
			// ADD-BY-LEETEN 12/30/2012-END
		}

		// ADD-BY-LEETEN 12/28/2012-BEGIN
		virtual	
		void
		_GetInteger(
			int eName,
			long *plValue,
			void* _Reserved = NULL
		)
		{
			CDecoderBase<ST, BT>::_GetInteger(eName, plValue);
			CSATSepDWTNetCDF::_GetInteger(eName, plValue);
			CSepDWTHeader::_GetInteger(eName, plValue);
		}

		// ADD-BY-LEETEN 12/28/2012-END

		virtual
		void
		_Allocate(
					void *_Reserved = NULL
					)
		{
			// ADD-BY-LEETEN 12/29/2012-BEGIN
			vusCachedNextOffsets.resize(this->uNrOfCoefs);
			// vusCachedBins.resize(this->uNrOfCoefs);
			// vCachedValues.resize(this->uNrOfCoefs);

			vuGlobalValueBase.resize(this->uNrOfUpdatingCoefs);
			vuLocalValueCount.resize(this->uNrOfUpdatingCoefs);
			size_t uValueBase = 0;
			for(size_t c = 0, w = 0; w < uNrOfUpdatingCoefs; w++)
			{
				size_t uNrOfLocalCoefs;
				vector<size_t> vuGlobalCoefBase;
				vector<size_t> vuLocalCoefLength;
				this->_ConvertWaveletToLevels(w, vuGlobalCoefBase, vuLocalCoefLength, uNrOfLocalCoefs);
				size_t uNrOfLocalValues = 0;
				for(size_t lc = 0; lc < uNrOfLocalCoefs; lc++, c++)
				{
					size_t gc = this->vuMapLocalToGlobal[c];
					uNrOfLocalValues += (size_t)this->vusCoefCounts[gc];
				}
				vuGlobalValueBase[w] = uValueBase;
				vuLocalValueCount[w] = uNrOfLocalValues;
				uValueBase += uNrOfLocalValues;
			}
			// ADD-BY-LEETEN 12/29/2012-END

			vCoefBins.resize(UGetNrOfBins());
			vCoefValues.resize(UGetNrOfBins());

			// sort the coef. indices by its # coefficients
			vector< pair<BT, long long> > vpairCoefCountIndex;
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
			// vector<size_t> vuCoefSub;	vuCoefSub.resize(UGetNrOfDims());
			vector<size_t> vuLevel;		vuLevel.resize(UGetNrOfDims());
			vector<size_t> vuGlobalBase;	vuGlobalBase.resize(UGetNrOfDims());
			vector<size_t> vuLocalLengths;	vuLocalLengths.resize(UGetNrOfDims());
			vector<size_t> vuLocalCoef;		vuLocalCoef.resize(UGetNrOfDims());
			vector<size_t> vuGlobalCoef;	vuGlobalCoef.resize(UGetNrOfDims());
			size_t uNrOfValuesInCore = 0;

			reverse(vpairCoefCountIndex.begin(), vpairCoefCountIndex.end());

			for(size_t c = 0; c < uNrOfCoefs && uNrOfValuesInCore <= uMaxNrOfValuesInCore; c++)
			{
				size_t uCoef = abs(vpairCoefCountIndex[c].second);
				if( vbFlagsCoefInCore[uCoef] )
					continue;

				_ConvertIndexToLevels
				(
					uCoef,
					vuLevel,
					vuLocalCoef,
					vuGlobalBase,
					vuLocalLengths
				);
				size_t uWavelet = UConvertSubToIndex(vuLevel, vuDimLevels);
				if( uNrOfValuesInCore + this->vuLocalValueCount[uWavelet] <= uMaxNrOfValuesInCore )
				{
					size_t uNrOfLocalCoefs;
					_ConvertWaveletToLevels(uWavelet, vuGlobalBase, vuLocalLengths, uNrOfLocalCoefs);

					for(size_t lc = 0; lc < uNrOfLocalCoefs; lc++)
					{
						_ConvertIndexToSub(lc, vuLocalCoef, vuLocalLengths);
						for(size_t d = 0; d < UGetNrOfDims(); d++)
							vuGlobalCoef[d] = vuGlobalBase[d] + vuLocalCoef[d];
						size_t uGlobalCoef = UConvertSubToIndex(vuGlobalCoef, vuCoefLengths);
						vbFlagsCoefInCore[uGlobalCoef] = true;
					}
					uNrOfValuesInCore += this->vuLocalValueCount[uWavelet];
				}
			}
			LOG_VAR(uNrOfNonZeroValues);
			LOG_VAR(uMaxNrOfValuesInCore);
			LOG_VAR(uNrOfValuesInCore);

			this->vpcCoefPools.resize(this->uNrOfUpdatingCoefs);	// allocate the pools

			size_t uNrOfQueries = (size_t) 1 << UGetNrOfDims();
			vpairLocalWaveletQueues.resize(uNrOfQueries * this->uNrOfUpdatingCoefs);
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
		  // ADD-BY-LEETEN 01/02/2013-BEGIN
		  const char* szFilepath = this->szFilepath;
		  // ADD-BY-LEETEN 01/02/2013-END

			// read the file size
			// ADD-BY-LEETEN 12/28/2012-BEGIN
			#if		WITH_BOOST
			size_t uFileSize = fs::file_size( szFilepath );
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

			// ADD-BY-LEETEN 01/05/2012-BEGIN
			if(uNrOfRangeQueries)
			{
				LOG_VAR(uMinNrOfCoefQueries);
				LOG_VAR(uMaxNrOfCoefQueries);
				LOG_VAR(uAccumNrOfCoefQueries/uNrOfRangeQueries);
			}
			// ADD-BY-LEETEN 01/05/2012-END
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
			// ADD-BY-LEETEN 01/27/2013-BEGIN
			LIBCLOCK_INIT(bIsPrintingTiming, __FUNCTION__);
			LIBCLOCK_BEGIN(bIsPrintingTiming);
			// ADD-BY-LEETEN 01/27/2013-END

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
			  for(size_t t = 0; t < (size_t)NR_OF_DIM_TYPES; t++)
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

			_Set(vuDimLengths, (BT)uNrOfBins);

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
			vector<TYPE_COEF_OFFSET>	vCoefOffsets;
			vector<TYPE_COEF_COUNT>		vCoefCounts;

			// now read the entire header
			vuCoefOffsets.resize(uNrOfCoefs);
			vusCoefCounts.resize(uNrOfCoefs);
			vCoefOffsets.resize(vuCoefOffsets.size());
			ASSERT_NETCDF(nc_get_var(
				iNcId,
				ncVarCoefOffset,
				vCoefOffsets.data()));
			for(size_t h = 0; h < vuCoefOffsets.size(); h++)
				vuCoefOffsets[h] = (size_t)vCoefOffsets[h];

			vCoefCounts.resize(vusCoefCounts.size());

			ASSERT_NETCDF(nc_get_var(
				iNcId,
				ncVarCoefCount,
				vCoefCounts.data()));
			for(size_t h = 0; h < vusCoefCounts.size(); h++)
				vusCoefCounts[h] = (BT)vCoefCounts[h];

			// ADD-BY-LEETEN 01/27/2013-BEGIN
			LIBCLOCK_END(bIsPrintingTiming);

			LIBCLOCK_BEGIN(bIsPrintingTiming);
			// ADD-BY-LEETEN 01/27/2013-END

			// ADD-BY-LEETEN 12/26/2012-BEGIN
			_Allocate();
			size_t uNrOfAllocatedPools = 0;
			// ADD-BY-LEETEN 12/26/2012-END

			// ADD-BY-LEETEN 01/27/2013-BEGIN
			LIBCLOCK_END(bIsPrintingTiming);

			LIBCLOCK_BEGIN(bIsPrintingTiming);
			// ADD-BY-LEETEN 01/27/2013-END
			/////////////////////////////////////////////////////////////////
			// now load the coefficients that can be in core
			vector<size_t> vuLocalCoefLengths;
			vuLocalCoefLengths.resize(UGetNrOfDims());
			vector<size_t> vuGlobalCoefBase;
			vuGlobalCoefBase.resize(UGetNrOfDims());

			vector<TYPE_COEF_BIN> vCoefBins;
			vCoefBins.resize(uNrOfNonZeroValues);
			ASSERT_NETCDF(nc_get_var(
				iNcId,
				ncVarCoefBin,
				vCoefBins.data()));

			vector<TYPE_COEF_VALUE> vCoefValues;
			vCoefValues.resize(uNrOfNonZeroValues);
			ASSERT_NETCDF(nc_get_var(
				iNcId,
				ncVarCoefValue,
				vCoefValues.data()));

			for(size_t i = 0, c = 0; c < this->uNrOfUpdatingCoefs; c++)
			{
				size_t uNrOfLocalCoefs; 
				this->_ConvertWaveletToLevels(c, vuGlobalCoefBase, vuLocalCoefLengths, uNrOfLocalCoefs);

				// ADD-BY-LEETEN 2013/08/11-BEGIN
				vector<size_t> vuPoolLevels;
				_ConvertIndexToSub(c, vuPoolLevels, vuDimLevels);

				// decide the wavelet weight
				double dWavelet = +1.0;

				vector<size_t> vuWaveletLengths;
				vuWaveletLengths.resize(UGetNrOfDims());

				for(size_t d = 0; d < vuWaveletLengths.size(); d++)
				{
					size_t uLevel = vuPoolLevels[d];
					if( uLevel >= 1 )
						dWavelet *= (double)(1 << (uLevel - 1));
					vuWaveletLengths[d] = (size_t)1<<(( !uLevel )?(vuDimLevels[d] - 1):(vuDimLevels[d] - uLevel));
				}
				dWavelet = sqrt(dWavelet);
				double dWeight = dWavelet/dWaveletDenomiator;
				// ADD-BY-LEETEN 2013/08/11-END


				// scan through all basis
				for(size_t valuei = 0, lc = 0; lc < uNrOfLocalCoefs; lc++, i++)
				{
					size_t uCoefIndex = this->vuMapLocalToGlobal[i];
					bool bIsInCore = this->vbFlagsCoefInCore[uCoefIndex];
					if( !vpcCoefPools[c] )
					{
						// MOD-BY-LEETEN 2013/10/27-FROM:						vpcCoefPools[c]= new CSepDWTPool<WT, BT>;
						vpcCoefPools[c]= new CSepDWTPool();
						// MOD-BY-LEETEN 2013/10/27-END
						vpcCoefPools[c]->_Set(
							(BT)UGetNrOfBins(),
							vuLocalCoefLengths,
							vuMaxCounts[c],
							true);

						// ADD-BY-LEETEN 2013/08/11-BEGIN
						vpcCoefPools[c]->_SetWaveletWeight(dWeight);
						vpcCoefPools[c]->_SetDataDimLengths(vuDimLengths);
						vpcCoefPools[c]->_SetWaveletLengths(vuWaveletLengths);
						// ADD-BY-LEETEN 2013/08/11-END
						uNrOfAllocatedPools++;
					}

					// ADD-BY-LEETEN 2013/08/11-BEGIN
					vector<size_t> vuSub;
					_ConvertIndexToSub(lc, vuSub, vuLocalCoefLengths);
					bool bIsOutOfDataDomain = false;
					for(size_t d = 0; d < vuSub.size(); d++) 
					{
						if( vuSub[d] >= (size_t)ceilf((float)vuDimLengths[d]/(float)vuWaveletLengths[d]) ) 
						{
							bIsOutOfDataDomain = true;
							break;
						}
					}

					if( bIsOutOfDataDomain )
						continue;
					// ADD-BY-LEETEN 2013/08/11-END
					this->vpcCoefPools[c]->_Copy(
						lc,
						vusCoefCounts[uCoefIndex], 
						&vCoefBins.data()[vuCoefOffsets[uCoefIndex]],
						&vCoefValues.data()[vuCoefOffsets[uCoefIndex]]);
				}
				if( vpcCoefPools[c] )	
					this->vpcCoefPools[c]->_Finalize(1.0);
			}

			// ADD-BY-LEETEN 12/26/2012-BEGIN
			LOG_VAR(uNrOfAllocatedPools);
			// ADD-BY-LEETEN 12/26/2012-END
			// #endif	// #if WITH_NETCDF 

			// ADD-BY-LEETEN 01/27/2013-BEGIN
			LIBCLOCK_END(bIsPrintingTiming);

			LIBCLOCK_PRINT(bIsPrintingTiming);
			// ADD-BY-LEETEN 01/27/2013-END
		}

		// ADD-BY-LEETEN 01/05/2013-BEGIN
		virtual
		void
		_GetCoefSparse
		(
			size_t uWavelet,
			size_t uLocalCoef,
			vector< pair<BT, WT> >& vpairCoefBinValues,
			void* _Reserved = NULL
		) 
		{
			// ADD-BY-LEETEN 2013/10/30-BEGIN
			if( uWavelet >= this->vpcCoefPools.size() ) {
				LOG_ERROR(fprintf(stderr, "Error: invalid index."));
				return;
			}
			// ADD-BY-LEETEN 2013/10/30-END

			if( this->vpcCoefPools[uWavelet] )
			{
				this->vpcCoefPools[uWavelet]->_GetCoefSparse
				(
					uLocalCoef,
					vpairCoefBinValues
				);
			}
			else
			{
				vector<size_t> vuLocalCoef;
				_ConvertIndexToSub(uLocalCoef, vuLocalCoef, this->vvuLocalLengths[uWavelet]);
				vector<size_t> vuGlobalCoef;
				vuGlobalCoef.resize(UGetNrOfDims());
				for(size_t d = 0; d < vuGlobalCoef.size(); d++)
					vuGlobalCoef[d] = this->vvuGlobalBase[uWavelet][d] + vuLocalCoef[d];
				size_t uGlobalCoefIndex = UConvertSubToIndex(vuGlobalCoef, this->vuCoefLengths);

				size_t uStart = (size_t)vuCoefOffsets[uGlobalCoefIndex];
				size_t uCount = (size_t)vusCoefCounts[uGlobalCoefIndex];

				vpairCoefBinValues.clear();
				if( uCount )
				{
					ASSERT_NETCDF(nc_get_vara(
						iNcId,
						ncVarCoefBin,
						&uStart,
						&uCount,
						vCoefBins.data()));
	
					ASSERT_NETCDF(nc_get_vara(
						iNcId,
						ncVarCoefValue,
						&uStart,
						&uCount,
						vCoefValues.data()));
					for(size_t i = 0; i < uCount; i++)
						vpairCoefBinValues.push_back(pair<BT, WT>((BT)vCoefBins[i], (WT)vCoefValues[i]));
				}	
				sort(vpairCoefBinValues.begin(), vpairCoefBinValues.end());			// ADD-BY-LEETEN 2013/11/02/2013
			}
		}
		// ADD-BY-LEETEN 01/05/2013-END

		// ADD-BY-LEETEN 2013/10/30-BEGIN
		virtual
		void
		_GetCoefSparse
		(
			const vector<size_t>& vuWaveletSub,
			const vector<size_t>& vuLocal,
			vector< pair<BT, WT> >& vpairCoefBinValues,
			void* _Reserved = NULL
		) 
		{
			vector<size_t> vuGlobalCoefBase;
			vector<size_t> vuLocalCoefLengths;
			_ConvertWaveletSubToLevels(vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths);

			_GetCoefSparse(
				UConvertSubToIndex(vuWaveletSub, vuDimLevels),
				UConvertSubToIndex(vuLocal, vuLocalCoefLengths), 
				vpairCoefBinValues);
		}

		virtual
		void
		_GetCoefSparse
		(
			size_t uWavelet,
			const vector<size_t>& vuLocal,
			vector< pair<BT, WT> >& vpairCoefBinValues,
			void* _Reserved = NULL
		) 
		{
			vector<size_t> vuGlobalCoefBase;
			vector<size_t> vuLocalCoefLengths;
			size_t uNrOfLocalCoefs;
			_ConvertWaveletToLevels(uWavelet, vuGlobalCoefBase, vuLocalCoefLengths, uNrOfLocalCoefs);
			size_t uLocal = UConvertSubToIndex(vuLocal, vuLocalCoefLengths);
			_GetCoefSparse(uWavelet, uLocal, vpairCoefBinValues);
		}

		virtual
		void
		_GetCoefSparse
		(
			const vector<size_t>& vuWaveletSub,
			size_t uLocal,
			vector< pair<BT, WT> >& vpairCoefBinValues,
			void* _Reserved = NULL
		) 
		{
			size_t uWavelet = UConvertSubToIndex(vuWaveletSub, vuDimLevels);
			_GetCoefSparse(uWavelet, uLocal, vpairCoefBinValues);
		}
		// ADD-BY-LEETEN 2013/10/30-END

		//! Return the sum of all bins at the given position
		virtual	
		void
		_GetAllSums
		(
			const vector<size_t>& vuPos,
			vector<ST>& vdSums,
			void *_Reserved = NULL
		)
		{
			vector<size_t> vuEmpty;	// this empty vector is used to speed up the access of elements
			vdSums.assign(UGetNrOfBins(), (ST)0);

			vector<size_t> vuLevels;			vuLevels.resize( UGetNrOfDims() );
			vector<size_t> vuGlobalCoefBase;	vuGlobalCoefBase.resize( UGetNrOfDims() );
			vector<size_t> vuLocalCoefSub;		vuLocalCoefSub.resize( UGetNrOfDims() );
			vector<size_t> vuCoefSub;			vuCoefSub.resize( UGetNrOfDims() );
		
			vector<double> vdWaveletBasis;
			vector<size_t> vuSubs;
			_GetBackwardWavelet(
				vuPos, 
				vuSubs, 
				vdWaveletBasis, 
				true);

			// ADD-BY-LEETEN 12/28/2012-BEGIN
			this->uNrOfQueries++;
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

				vector< pair<BT, WT> > vpairCoefBinValues;
				size_t uLocalCoef = UConvertSubToIndex(vuLocalCoefSub, this->vvuLocalLengths[c]);
				this->_GetCoefSparse
				(
					c,
					uLocalCoef,
					vpairCoefBinValues
				);
				for(typename vector< pair<BT, WT> >::iterator
					ivpairCoefs = vpairCoefBinValues.begin();
					ivpairCoefs != vpairCoefBinValues.end();
					ivpairCoefs++ )
					vdSums[ivpairCoefs->first] += ivpairCoefs->second * dWavelet;
			}
			// ADD-BY-LEETEN 12/28/2012-BEGIN
			this->uAccumNrOfIORequest += uNrOfIORequest;
			this->uMaxNrOfIORequest = max(this->uMaxNrOfIORequest, uNrOfIORequest);
			if( uNrOfIORequest )
			{
				if( !this->uMinNrOfIORequest )
					this->uMinNrOfIORequest = uNrOfIORequest;
				else
					this->uMinNrOfIORequest = min(this->uMinNrOfIORequest, uNrOfIORequest);
			}
			// ADD-BY-LEETEN 12/28/2012-END
			for(size_t b = 0; b < UGetNrOfBins(); b++)
				vdSums[b] /= dWaveletDenomiator;
		}

		// ADD-BY-LEETEN 01/05/2013-BEGIN
		//! Return the sum of all bins at the given position
		virtual	
		void
		_GetRegionSums
		(
			const vector<size_t>& vuLeft,
			const vector<size_t>& vuRight,
			vector<ST>& vdSums,
			void *_Reserved = NULL
		)
		{
			vdSums.assign(UGetNrOfBins(), (ST)0);

			size_t  uNrOfCoefQueries = 0;
			#if	!MERGE_COEFFICIENTS_PER_REGION
			size_t uNrOfQueries = (size_t)1 << this->UGetNrOfDims();
			vector<size_t> vuQueryPos;
			vuQueryPos.resize(this->UGetNrOfDims());
			vector<ST> vQuerySums;
			for(size_t q = 0; q < uNrOfQueries; q++)
			{
				int iSign = 1;
				for(size_t d = 0, j = q; d < this->UGetNrOfDims(); d++, j /= 2)
				{
					vuQueryPos[d] = (j % 2)?vuLeft[d]:vuRight[d];
					iSign *= (j % 2)?(-1):(+1);
				}
				_GetAllSums(vuQueryPos, vQuerySums);
				for(size_t b = 0; b < vdSums.size(); b++)
					vdSums[b] += (ST)iSign * vQuerySums[b];
			}
			uNrOfCoefQueries = uNrOfQueries * this->uNrOfUpdatingCoefs;
			#else	// #if	!MERGE_COEFFICIENTS_PER_REGION
			size_t uNrOfQueries = (size_t)1 << this->UGetNrOfDims();
			vector<size_t> vuQueryPos;
			vuQueryPos.resize(this->UGetNrOfDims());
			vector<size_t> vuLocalCoefSub;
			vuLocalCoefSub.resize(this->UGetNrOfDims());
			for(size_t q = 0; q < uNrOfQueries; q++)
			{
				int iSign = 1;
				for(size_t d = 0, j = q; d < this->UGetNrOfDims(); d++, j /= 2)
				{
					vuQueryPos[d] = (j % 2)?vuLeft[d]:vuRight[d];
					iSign *= (j % 2)?(-1):(+1);
				}

				// find the corresponding wavelet coefficient
				vector<double> vdWaveletBasis;
				vector<size_t> vuSubs;
				_GetBackwardWavelet(
					vuQueryPos, 
					vuSubs, 
					vdWaveletBasis, 
					true);

				// now find the combination of the coefficients of all dimensions 
				for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
				{
					double dWavelet = 1.0;
					for(size_t d = 0, uBase = 0;
						d < UGetNrOfDims(); 
						uBase += vuDimMaxLevels[d], d++, p++)
					{
						// compute the wavelet
						dWavelet *= vdWaveletBasis[uBase + vuCoefDim2Level[p]];	

						// compute the subscript and size for the current basis
						vuLocalCoefSub[d] = vuSubs[uBase + vuCoefDim2Level[p]];
					}
					size_t uLocalCoef = UConvertSubToIndex(vuLocalCoefSub, vvuLocalLengths[c]);

					this->vpairLocalWaveletQueues[c * uNrOfQueries + q] = pair<size_t, WT>(uLocalCoef, (WT)iSign * (WT)dWavelet);
				}
			}

			for(size_t i = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				sort(	vpairLocalWaveletQueues.begin() + i,
						vpairLocalWaveletQueues.begin() + i + uNrOfQueries);
				pair<size_t, WT> vpairPrevLocalCoefQueues = pair<size_t, WT>(vpairLocalWaveletQueues[i].first, (WT)0);
				for(size_t q = 0; q < uNrOfQueries + 1; i += (q < uNrOfQueries)?1:0, q++)
				{
					bool bIsNew = false;
					if( q == uNrOfQueries )
						bIsNew = true;
					else 
					{
						// LOG(printf("q(%d): %d %f", q, vpairLocalWaveletQueues[i].first, vpairLocalWaveletQueues[i].second));
						if( vpairPrevLocalCoefQueues.first	== vpairLocalWaveletQueues[i].first )
						{
							vpairPrevLocalCoefQueues.second += vpairLocalWaveletQueues[i].second;
							// LOG_VAR(vpairPrevLocalCoefQueues.second);
						}
						else
							bIsNew = true;
					}

					if( bIsNew )
					{
						if(	this->dWaveletThreshold < fabs((double)vpairPrevLocalCoefQueues.second) )
						{
							vector< pair<BT, WT> > vpairCoefBinValues;
							this->_GetCoefSparse
							(
								c,
								vpairPrevLocalCoefQueues.first,
								vpairCoefBinValues
							);
							for(typename vector< pair<BT, WT> >::iterator
								ivpairCoefs = vpairCoefBinValues.begin();
								ivpairCoefs != vpairCoefBinValues.end();
								ivpairCoefs++ )
								vdSums[ivpairCoefs->first] += ivpairCoefs->second * vpairPrevLocalCoefQueues.second;
							uNrOfCoefQueries++;
						}
						if( i < vpairLocalWaveletQueues.size() )	// ADD-BY-LEETEN 05/05/2013
						vpairPrevLocalCoefQueues = vpairLocalWaveletQueues[i];
					}
				}
			}
			for(size_t b = 0; b < UGetNrOfBins(); b++)
				vdSums[b] /= dWaveletDenomiator;
			#endif	// #if	!MERGE_COEFFICIENTS_PER_REGION
			// LOG_VAR();
			// ADD-BY-LEETEN 01/05/2012-BEGIN
			if( uNrOfCoefQueries )
				uMinNrOfCoefQueries = ( !uMinNrOfCoefQueries )?uNrOfCoefQueries:min(uMinNrOfCoefQueries, uNrOfCoefQueries);
			uMaxNrOfCoefQueries = max(uMaxNrOfCoefQueries, uNrOfCoefQueries);
			uAccumNrOfCoefQueries += uNrOfCoefQueries;
			uNrOfRangeQueries++;
			// ADD-BY-LEETEN 01/05/2012-END
		}
		// ADD-BY-LEETEN 01/05/2013-END

		// ADD-BY-LEETEN 01/18/2012-BEGIN
		virtual	
		void
		_GetDecodedSize
		(
			vector<size_t>& vuDecodedSize,
			void *_Reserved = NULL
		) const
		{
			if(vuDecodedSize.size() != UGetNrOfDims())
				vuDecodedSize.resize(UGetNrOfDims());
			copy(vuCoefLengths.begin(), vuCoefLengths.end(), vuDecodedSize.begin());
		}
		// ADD-BY-LEETEN 01/18/2012-END

		// ADD-BY-LEETEN 12/29/2012-BEGIN
		virtual
		void
		_DecodeBin
		(
			const BT& usBin,
			vector<ST>& vSAT,
			void *_Reserved = NULL
		)
		{
		  // ADD-BY-LEETEN 01/02/2013-BEGIN
		  const bool bIsPrintingDecodeBinTiming = this->bIsPrintingDecodeBinTiming;
		  // ADD-BY-LEETEN 01/02/2013-END
			if( uNrOfCoefs != vSAT.size() )
				vSAT.resize(uNrOfCoefs);

			LIBCLOCK_INIT(	bIsPrintingDecodeBinTiming, __FUNCTION__);
			LIBCLOCK_BEGIN(	bIsPrintingDecodeBinTiming);
			vector<size_t> vuEmpty;
			vector<size_t> vuGlobalCoefBase, vuLocalCoefLengths;

			for(size_t c = 0, w = 0; w < uNrOfUpdatingCoefs; w++)
			{
				vector<TYPE_COEF_BIN>	vCoefBins;		
				vector<TYPE_COEF_VALUE> vCoefValues;	
				if( !this->vpcCoefPools[w] )
				{
					size_t uLocalValueBase = this->vuGlobalValueBase[w];
					size_t uNrOfLocalValues = this->vuLocalValueCount[w];
					vCoefBins.resize(uNrOfLocalValues);
					vCoefValues.resize(uNrOfLocalValues);

					ASSERT_NETCDF(nc_get_vara(
						iNcId,
						ncVarCoefBin,
						&uLocalValueBase,
						&uNrOfLocalValues,
						vCoefBins.data()));
	
					ASSERT_NETCDF(nc_get_vara(
						iNcId,
						ncVarCoefValue,
						&uLocalValueBase,
						&uNrOfLocalValues,
						vCoefValues.data()));
				}
				///////////////////////////////////////
				size_t uNrOfLocalCoefs = 1;
				_ConvertWaveletToLevels(w, vuGlobalCoefBase, vuLocalCoefLengths, uNrOfLocalCoefs);
				for(size_t lc = 0,		uLocalValueBase = 0, gc = 0;
					lc < uNrOfLocalCoefs; 
					lc++,				uLocalValueBase += (size_t)vusCoefCounts[gc], c++)
				{
					gc = vuMapLocalToGlobal[c];

					BT usCount = vusCoefCounts[gc];
					BT usNextOffset = vusCachedNextOffsets[gc];
					BT usFetchedBin; // = vusCachedBins[gc];
					WT FetchedValue; // = vCachedValues[gc];

					// move the cache bin till it is great than or equal to the given bin
					if( !this->vbFlagsCoefInCore[gc] )
					{
						for(; usNextOffset < usCount; usNextOffset++)
						{
							// if( !usNextOffset || usBin < usFetchedBin)  
							usFetchedBin = vCoefBins[uLocalValueBase + (size_t)usNextOffset];
							FetchedValue = (WT)vCoefValues[uLocalValueBase + (size_t)usNextOffset];
							if( usFetchedBin >= usBin )
								break;
						}
					}
					else
					{
						for(; usNextOffset < usCount; usNextOffset++)
						{
							// if( !usNextOffset || usBin < usFetchedBin) 
							this->vpcCoefPools[w]->_GetAtOffset(usNextOffset, vuEmpty, lc, usFetchedBin, FetchedValue);
							if( usFetchedBin >= usBin )
								break;
						}
					}

					// if the cached bin is equal to the given bin, the value is the cached one
					vSAT[gc] = ( 0 < usCount && usNextOffset < usCount && usFetchedBin == usBin )?FetchedValue:(WT)0;
					// update the cache 
					vusCachedNextOffsets[gc] = usNextOffset;
				}
			}
			LIBCLOCK_END(bIsPrintingDecodeBinTiming);

			//////////////////////////////////////////////////
			// now apply IDWT
			LIBCLOCK_BEGIN(bIsPrintingDecodeBinTiming);
		      for(size_t uOffset = 1, d = 0, uCoefLength = 1; 
				d < UGetNrOfDims(); 
				uOffset *= uCoefLength, d++)
			{
			  uCoefLength = this->vuCoefLengths[d]; // ADD-BY-LEETEN 12/31/2012

				if( 1 == uCoefLength )
					continue;

				vector<ST> pvDsts[2];
				for(size_t i = 0; i < 2; i++)
					pvDsts[i].resize(uCoefLength);

				/*
				vector<size_t> vuScanLineIndices;
				vuScanLineIndices.resize(uCoefLength);
				*/

				size_t uNrOfLevels = vuDimLevels[d] - 1;
				size_t uNrOfScanLines = this->uNrOfCoefs / uCoefLength;

				for(size_t i = 0; i < uNrOfScanLines; i++)
				{
					size_t uScanlineBase = vvuSliceScanlineBase[d][i];
					ST *vSrc = &vSAT.data()[uScanlineBase];
					pvDsts[0][0] = vSrc[0];
					_IDWT1D(
						vSrc, 
						uOffset,
						pvDsts[0].data(), 
						pvDsts[1].data(), 
						vuCoefLengths[d],
						vuDimLengths[d],
						2, 
						uNrOfLevels - 1);
					vector<ST>& vDST = pvDsts[uNrOfLevels%2];
					for(size_t si = 0, di = 0; di < uCoefLength; di++, si += uOffset)
							vSrc[si] = vDST[di];
				}
			}
			LIBCLOCK_END(bIsPrintingDecodeBinTiming);
			LIBCLOCK_PRINT(bIsPrintingDecodeBinTiming);
		}
		// ADD-BY-LEETEN 12/29/2012-END

		// ADD-BY-LEETEN 01/02/2013-BEGIN
		virtual
		void
		_ClampToDataSize(
			const vector<ST>& vCoefField,
			vector<ST>& vDataField,
			void* _Reserved = NULL
			)
		{
			// only keep the entropy field within the data range
			if( uDataSize != vDataField.size() )
				vDataField.resize(uDataSize);
			vector<size_t> vuSub;
			for(size_t d = 0; d < uDataSize; d++)
			{
				_ConvertIndexToSub(d, vuSub, vuDimLengths);
				vDataField[d] = vCoefField[UConvertSubToIndex(vuSub, vuCoefLengths)];
			}
		}

		virtual
		void
		_ClampBorder(
			vector<ST>& vField,	
			const vector<int>& viLeft, 
			const vector<int>& viRight, 
			void* _Reserved = NULL
			)
		{
			vector<size_t> vuSub;
			for(size_t d = 0; d < uDataSize; d++)
			{
				vector<size_t> vuSub;
				_ConvertIndexToSub(d, vuSub, vuDimLengths);
				bool bIsNearBorder = false;
				for(size_t dim = 0; dim < this->UGetNrOfDims(); dim++)
					if( 0 > (int)vuSub[dim] + viLeft[dim] || 
							(int)vuSub[dim] + viLeft[dim] >= vuDimLengths[dim] ||
						0 > (int)vuSub[dim] + viRight[dim] || 
							(int)vuSub[dim] + viRight[dim] >= vuDimLengths[dim] )
					{
						bIsNearBorder = true;
						break;
					}

				if( bIsNearBorder )
				{
					vField[d] = (ST)0;
					continue;
				}
			}
		}
		// ADD-BY-LEETEN 01/02/2013-END

		// ADD-BY-LEETEN 01/23/2013-BEGIN
		//! Apply statistics to the wavelet coefficients in the specified wavelet function
		virtual
		void
		_CompBlockStatistics
		(
			size_t uBlockLevel,	//!< The subscvript of the wavelet 
			ST (*Stat)(const vector< pair< BT, WT > >&),
			const char *szStaticsFilepath,
			void* _Reserved = NULL
		)
		{		
			if( 0 == uBlockLevel )
				return;

			vector<size_t> vuWaveletSub;
			vuWaveletSub.resize(UGetNrOfDims());
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				if( vuDimLevels[d] < uBlockLevel )
					return;
				else
					vuWaveletSub[d] = vuDimLevels[d] - uBlockLevel;

			LOG_VAR(szStaticsFilepath);

			vector<size_t> vuGlobalCoefBase;
			vector<size_t> vuLocalCoefLengths;
			this->_ConvertWaveletSubToLevels(
				vuWaveletSub,
				vuGlobalCoefBase,
				vuLocalCoefLengths);
			size_t uNrOfLocalCoefs = 1;
			for(size_t d = 0; d < this->UGetNrOfDims(); d++)
				uNrOfLocalCoefs *= vuLocalCoefLengths[d];

			// now skip the portion outside the data
			size_t uLocalDataLength = 1;
			vector<size_t> vuLocalDataLengths;
			vuLocalDataLengths.resize(UGetNrOfDims());
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				vuLocalDataLengths[d] = (size_t)ceil( (double)vuDimLengths[d] / (double)((size_t)1<<uBlockLevel) );
				uLocalDataLength *= vuLocalDataLengths[d];
			}

			// scan through all coefficients to compute the statistics
			vector<ST> vStat;
			vStat.resize(uLocalDataLength);

			size_t uWavelet = UConvertSubToIndex(vuWaveletSub, vuDimLevels);

			vector< pair< BT, WT > > vpairCoefBinValues;
			vector<size_t> vuLocalDataSub;
			for(size_t l = 0; l < uLocalDataLength; l++)
			{
				_ConvertIndexToSub(l, vuLocalDataSub, vuLocalDataLengths);
				this->_GetCoefSparse(uWavelet, UConvertSubToIndex(vuLocalDataSub, vuLocalCoefLengths), vpairCoefBinValues);
				vStat[l] = Stat(vpairCoefBinValues);
			}

			// now save the file to NRRD format
			_SaveNrrd<ST>(vuLocalDataLengths, vStat.data(), szStaticsFilepath);
		}
		// ADD-BY-LEETEN 01/23/2013-END

		CSATSepDWTDecoder():
			// ADD-BY-LEETEN 01/05/2012-BEGIN
			uMinNrOfCoefQueries(0),
			uMaxNrOfCoefQueries(0),
			uAccumNrOfCoefQueries(0),
			uNrOfRangeQueries(0),
			// ADD-BY-LEETEN 01/05/2012-END
			CDecoderBase<ST, BT>(),
			CSATSepDWTNetCDF(),
			CSepDWTHeader()
		{
		}

		virtual	// ADD-BY-LEETEN 01/02/2013
		~CSATSepDWTDecoder()
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
