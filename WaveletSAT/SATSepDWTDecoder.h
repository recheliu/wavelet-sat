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
#include "DecoderBase.h"	// ADD-BY-LEETEN 01/02/2013

#include "liblog.h"	
#if	!WITH_SMART_PTR	// ADD-BY-LEETEN 12/30/2012
#include "libbuf.h"
// ADD-BY-LEETEN 12/30/2012-BEGIN
#else	// #if	!WITH_SMART_PTR
#include <boost/shared_array.hpp>
#endif	// #if	!WITH_SMART_PTR
// ADD-BY-LEETEN 12/30/2012-END

#if	WITH_NETCDF
#include <netcdf.h>
#include "lognc.h"
#endif	// #if	WITH_NETCDF

// ADD-BY-LEETEN 12/28/2012-BEGIN
// if this is non-0, the header with most coefficients will be place in core first
#define	IS_SELECTING_LONGEST_FIRST					1

// if this is non-0, when a coefficient is chosen, other coefficients of the same wavelet will be chosen as well
#define IS_SELECTING_THE_SAME_WAVELET			1
// ADD-BY-LEETEN 12/28/2012-END

// ADD-BY-LEETEN 01/05/2013-BEGIN
// If this is non-0, the coefficient for a region query will be merged to gether to reduce the access of coefficients
#define MERGE_COEFFICIENTS_PER_REGION			1
// ADD-BY-LEETEN 01/05/2013-END

// ADD-BY-LEETEN 01/05/2013-BEGIN
#if		MERGE_COEFFICIENTS_PER_REGION	
// If this is non-0, the queue to hold the coefficient will be pre-allocated, whose size is (log n)^d * 2^d.
#define WITH_PRE_ALLOCATED_QUEUES				1
#endif	// #if		MERGE_COEFFICIENTS_PER_REGION	
// ADD-BY-LEETEN 01/05/2013-END

namespace WaveletSAT
{
	#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
	template<
		typename DT,	//!< Type of the data
		typename WT		//!< Type of the wavelet coefficients
	>
	#else	// MOD-BY-LEETEN 01/03/2013-TO:
	template<
		// typename DT,				//!< Type of the data
		typename ST = typeSum,		//!< Type of the sum
		typename BT = typeBin,		//!< Type of the bin
		typename WT = typeWavelet	//!< Type of the wavelet coefficientsd
	>
	#endif	// MOD-BY-LEETEN 01/03/2013-END
	// The class that access the coefficients from files (in NetCDF format)
	class CSATSepDWTDecoder:
		#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
		virtual public CDecoderBase<DT>,	// ADD-BY-LEETEN 01/02/2013
		virtual public CSATSepDWTNetCDF,
		virtual public CSepDWTHeader
		#else	// MOD-BY-LEETEN 01/03/2013-TO:
		virtual public CDecoderBase<ST, BT>,
		virtual public CSATSepDWTNetCDF,
		// MOD-BY-LEETEN 01/04/2013-FROM:		virtual public CSepDWTHeader<ST, BT>
		virtual public CSepDWTHeader
		// MOD-BY-LEETEN 01/04/2013-END
		#endif	// MOD-BY-LEETEN 01/03/2013-END
	{
protected:	
			#if	0	// DEL-BY-LEETEN 01/02/2013-BEGIN
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

			// ADD-BY-LEETEN 12/29/2012-BEGIN
			bool bIsPrintingDecodeBinTiming;
			#endif	// DEL-BY-LEETEN 01/02/2013-END

			// ADD-BY-LEETEN 01/05/2013-BEGIN
			#if		WITH_PRE_ALLOCATED_QUEUES	
			vector< pair<size_t, WT> > vpairLocalWaveletQueues; 
			#endif	// #if	WITH_PRE_ALLOCATED_QUEUES
			// ADD-BY-LEETEN 01/05/2013-END

			// ADD-BY-LEETEN 01/05/2012-BEGIN
			size_t uMinNrOfCoefQueries;
			size_t uMaxNrOfCoefQueries;
			size_t uAccumNrOfCoefQueries;
			size_t uNrOfRangeQueries;
			// ADD-BY-LEETEN 01/05/2012-END

			// MOD-BY-LEETEN 01/03/2013-FROM:			vector<unsigned short>	vusCachedNextOffsets;
			vector<BT>	vusCachedNextOffsets;
			// MOD-BY-LEETEN 01/03/2013-END

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
			// MOD-BY-LEETEN 01/03/2013-FROM:			vector<unsigned short>	vusCoefCounts;
			vector<BT>	vusCoefCounts;
			// MOD-BY-LEETEN 01/03/2013-END

			//! #Coef. values stored in full arrays
			size_t uMaxNrOfValuesInCore;	

			//! vector of pointers to the coefficient pools. 
			/*!
			If the poiner is NULL, it means that the coefficients are out of core.
			*/
			// MOD-BY-LEETEN 01/03/2013-FROM:			vector< CSepDWTPool<WT, unsigned short>* > vpcCoefPools;
			vector< CSepDWTPool<WT, BT>* > vpcCoefPools;
			// MOD-BY-LEETEN 01/03/2013-END

			#if	!WITH_SMART_PTR		// ADD-BY-LEETEN 12/30/2012
			//! An array to store the coefficient bins.
			TBuffer<TYPE_COEF_BIN>		pCoefBins;

			//! An array to store the coefficient values.
			TBuffer<TYPE_COEF_VALUE>	pCoefValues;
			// ADD-BY-LEETEN 12/30/2012-BEGIN
			#else	// #if	!WITH_SMART_PTR	
			boost::shared_array<TYPE_COEF_BIN>		pCoefBins;

			boost::shared_array<TYPE_COEF_VALUE>	pCoefValues;
			#endif	// #if	!WITH_SMART_PTR	
			// ADD-BY-LEETEN 12/30/2012-END
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
			#if	0	// DEL-BY-LEETEN 01/02/2013-BEGIN
			ACCUM_NR_OF_IO_REQUESTS,
			MAX_NR_OF_IO_REQUESTS,
			MIN_NR_OF_IO_REQUESTS,
			RESET_IO_COUNTERS,
			
			//! Parameter specifying whether to print the timing for the method _Decodebin().
			PRINT_DECODE_BIN_TIMING,
			#endif	// DEL-BY-LEETEN 01/02/2013-END
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
			// MOD-BY-LEETEN 01/03/2013-FROM:			CDecoderBase<DT>::_SetInteger(eName, lValue);	// ADD-BY-LEETEN 01/02/2013
			CDecoderBase<ST, BT>::_SetInteger(eName, lValue);	// ADD-BY-LEETEN 01/02/2013
			// MOD-BY-LEETEN 01/03/2013-END
			CSATSepDWTNetCDF::_SetInteger(eName, lValue);
			// MOD-BY-LEETEN 01/03/2013-FROM:			CSepDWTHeader::_SetInteger(eName, lValue);
			// MOD-BY-LEETEN 01/04/2013-FROM:			CSepDWTHeader<ST, BT>::_SetInteger(eName, lValue);
			CSepDWTHeader::_SetInteger(eName, lValue);
			// MOD-BY-LEETEN 01/04/2013-END
			// MOD-BY-LEETEN 01/03/2013-END
			// ADD-BY-LEETEN 12/30/2012-BEGIN
			switch(eName)
			{
			#if 0 			  // MOD-BY-LEETEN 01/03/2013-FROM:
			case RESET_IO_COUNTERS:
				#if	0	// DEL-BY-LEETEN 01/02/2013-BEGIN
				uAccumNrOfIORequest = 0;
				uMaxNrOfIORequest = 0;
				#endif	// DEL-BY-LEETEN 01/02/2013-END
				uMinNrOfIORequest = uNrOfUpdatingCoefs;
				// DEL-BY-LEETEN 01/02/2013:	uNrOfQueries = 0;
				break;
			#else // MOD-BY-LEETEN 01/02/2013-TO:
			#if	0	// DEL-BY-LEETEN 01/05/2012-BEGIN
			case CDecoderBase<ST, BT>::RESET_IO_COUNTERS:
				this->uMinNrOfIORequest = uNrOfUpdatingCoefs;
				break;
			#endif	// DEL-BY-LEETEN 01/05/2012-END
			#endif // MOD-BY-LEETEN 01/02/2013-END
			#if	0	// DEL-BY-LEETEN 01/02/2013-BEGIN
			case PRINT_DECODE_BIN_TIMING:
				bIsPrintingDecodeBinTiming = (!lValue)?false:true;
				break;
			#endif	// DEL-BY-LEETEN 01/02/2013-END
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
			#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
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
			#else	// MOD-BY-LEETEN 01/02/2013-TO:
			// MOD-BY-LEETEN 01/03/2013-FROM:			CDecoderBase<DT>::_GetInteger(eName, plValue);
			CDecoderBase<ST, BT>::_GetInteger(eName, plValue);
			// MOD-BY-LEETEN 01/03/2013-END
			#endif	// MOD-BY-LEETEN 01/02/2013-END
			// ADD-BY-LEETEN 12/30/2012-BEGIN
			CSATSepDWTNetCDF::_GetInteger(eName, plValue);
			// MOD-BY-LEETEN 01/03/2013-FROM:			CSepDWTHeader::_GetInteger(eName, plValue);
			// MOD-BY-LEETEN 01/04/2013-FROM:			CSepDWTHeader<ST, BT>::_GetInteger(eName, plValue);
			CSepDWTHeader::_GetInteger(eName, plValue);
			// MOD-BY-LEETEN 01/04/2013-END
			// MOD-BY-LEETEN 01/03/2013-END
			// ADD-BY-LEETEN 12/30/2012-END
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
					// MOD-BY-LEETEN 01/03/2013-FROM:					uNrOfLocalValues += this->vusCoefCounts[gc];
					uNrOfLocalValues += (size_t)this->vusCoefCounts[gc];
					// MOD-BY-LEETEN 01/03/2013-END
				}
				vuGlobalValueBase[w] = uValueBase;
				vuLocalValueCount[w] = uNrOfLocalValues;
				uValueBase += uNrOfLocalValues;
			}
			// ADD-BY-LEETEN 12/29/2012-END

			#if	!WITH_SMART_PTR	// ADD-BY-LEETEN 12/30/2012
			pCoefBins.alloc(UGetNrOfBins());
			pCoefValues.alloc(UGetNrOfBins());
			// ADD-BY-LEETEN 12/30/2012-BEGIN
			#else	// #if	!WITH_SMART_PTR
			pCoefBins.reset(new TYPE_COEF_BIN[UGetNrOfBins()]);
			pCoefValues.reset(new TYPE_COEF_VALUE[UGetNrOfBins()]);
			#endif	// #if	!WITH_SMART_PTR
			// ADD-BY-LEETEN 12/30/2012-END
			// sort the coef. indices by its # coefficients
			// MOD-BY-LEETEN 01/03/2013-FROM:			vector< pair<unsigned short, long long> > vpairCoefCountIndex;
			vector< pair<BT, long long> > vpairCoefCountIndex;
			// MOD-BY-LEETEN 01/03/2013-END
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
				#endif	// #if		!IS_SELECTING_THE_SAME_WAVELET
			}
			LOG_VAR(uNrOfNonZeroValues);
			LOG_VAR(uMaxNrOfValuesInCore);
			LOG_VAR(uNrOfValuesInCore);

			this->vpcCoefPools.resize(this->uNrOfUpdatingCoefs);	// allocate the pools

			// ADD-BY-LEETEN 01/05/2013-BEGIN
			#if		WITH_PRE_ALLOCATED_QUEUES
			size_t uNrOfQueries = (size_t) 1 << UGetNrOfDims();
			vpairLocalWaveletQueues.resize(uNrOfQueries * this->uNrOfUpdatingCoefs);
			#endif	// #if	WITH_PRE_ALLOCATED_QUEUES
			// ADD-BY-LEETEN 01/05/2013-END
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

			// MOD-BY-LEETEN 01/03/2013-FROM:			_Set(vuDimLengths, uNrOfBins);
			_Set(vuDimLengths, (BT)uNrOfBins);
			// MOD-BY-LEETEN 01/03/2013-END

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
			#if	!WITH_SMART_PTR		// ADD-BY-LEETEN 12/30/2012
			TBuffer<TYPE_COEF_OFFSET> pCoefOffsets;
			TBuffer<TYPE_COEF_COUNT> pCoefCounts;
			#else	// #if	!WITH_SMART_PTR
			boost::shared_array<TYPE_COEF_OFFSET> pCoefOffsets;
			boost::shared_array<TYPE_COEF_COUNT> pCoefCounts;
			#endif	// #if	!WITH_SMART_PTR

			// now read the entire header
			vuCoefOffsets.resize(uNrOfCoefs);
			vusCoefCounts.resize(uNrOfCoefs);
			#if	!WITH_SMART_PTR		// ADD-BY-LEETEN 12/30/2012
			pCoefOffsets.alloc(vuCoefOffsets.size());
			// ADD-BY-LEETEN 12/30/2012-BEGIN
			#else	// #if	!WITH_SMART_PTR
			pCoefOffsets.reset(new TYPE_COEF_OFFSET[vuCoefOffsets.size()]);
			#endif	// #if	!WITH_SMART_PTR
			// ADD-BY-LEETEN 12/30/2012-END
			ASSERT_NETCDF(nc_get_var(
				iNcId,
				ncVarCoefOffset,
				(void*)&pCoefOffsets[0]));
			for(size_t h = 0; h < vuCoefOffsets.size(); h++)
				vuCoefOffsets[h] = (size_t)pCoefOffsets[h];
			#if	!WITH_SMART_PTR	// ADD-BY-LEETEN 12/30/2012
			pCoefOffsets.free();
			// ADD-BY-LEETEN 12/30/2012-BEGIN
			#else	// #if	!WITH_SMART_PTR
			pCoefOffsets.reset();
			#endif	// #if	!WITH_SMART_PTR
			// ADD-BY-LEETEN 12/30/2012-END

			#if	!WITH_SMART_PTR	// ADD-BY-LEETEN 12/30/2012
			pCoefCounts.alloc(vusCoefCounts.size());
			// ADD-BY-LEETEN 12/30/2012-BEGIN
			#else	// #if	!WITH_SMART_PTR
			pCoefCounts.reset(new TYPE_COEF_COUNT[vusCoefCounts.size()]);
			#endif	// #if	!WITH_SMART_PTR
			// ADD-BY-LEETEN 12/30/2012-END
			ASSERT_NETCDF(nc_get_var(
				iNcId,
				ncVarCoefCount,
				(void*)&pCoefCounts[0]));
			for(size_t h = 0; h < vusCoefCounts.size(); h++)
				// MOD-BY-LEETEN 01/03/2013-FROM:				vusCoefCounts[h] = (unsigned short)pCoefCounts[h];
				vusCoefCounts[h] = (BT)pCoefCounts[h];
				// MOD-BY-LEETEN 01/03/2013-END
			#if	!WITH_SMART_PTR	// ADD-BY-LEETEN 12/30/2012
			pCoefCounts.free();
			// ADD-BY-LEETEN 12/30/2012-BEGIN
			#else	// #if	!WITH_SMART_PTR
			pCoefCounts.reset();
			#endif	// #if	!WITH_SMART_PTR
			// ADD-BY-LEETEN 12/30/2012-END

			// ADD-BY-LEETEN 12/26/2012-BEGIN
			_Allocate();
			size_t uNrOfAllocatedPools = 0;
			// ADD-BY-LEETEN 12/26/2012-END

			/////////////////////////////////////////////////////////////////
			// now load the coefficients that can be in core
			vector<size_t> vuLocalCoefLengths;
			vuLocalCoefLengths.resize(UGetNrOfDims());
			vector<size_t> vuGlobalCoefBase;
			vuGlobalCoefBase.resize(UGetNrOfDims());

			size_t puStart[NC_MAX_DIMS];
			size_t puCount[NC_MAX_DIMS];
			for(size_t i = 0, c = 0; c < this->uNrOfUpdatingCoefs; c++)
			{
				size_t uNrOfLocalCoefs; 
				this->_ConvertWaveletToLevels(c, vuGlobalCoefBase, vuLocalCoefLengths, uNrOfLocalCoefs);
				for(size_t d = 0; d < UGetNrOfDims(); d++)
				{
					puStart[UGetNrOfDims() - 1 - d] = vuGlobalCoefBase[d];
					puCount[UGetNrOfDims() - 1 - d] = vuLocalCoefLengths[d];
				}

				// ADD-BY-LEETEN 12/29/2012-BEGIN
				#if	!WITH_SMART_PTR	// ADD-BY-LEETEN 12/30/2012
				TBuffer<TYPE_COEF_OFFSET> pLocalCoefOffsets;
				TBuffer<TYPE_COEF_COUNT> pLocalCoefCounts;
				// ADD-BY-LEETEN 12/29/2012-END
				pLocalCoefOffsets.alloc(uNrOfLocalCoefs);
				pLocalCoefCounts.alloc(uNrOfLocalCoefs);
				// ADD-BY-LEETEN 12/30/2012-BEGIN
				#else	// #if	!WITH_SMART_PTR	
				boost::shared_array<TYPE_COEF_OFFSET>	pLocalCoefOffsets(new TYPE_COEF_OFFSET[uNrOfLocalCoefs]);
				boost::shared_array<TYPE_COEF_COUNT>	pLocalCoefCounts(new TYPE_COEF_COUNT[uNrOfLocalCoefs]);;
				#endif	//	#if	!WITH_SMART_PTR	
				// ADD-BY-LEETEN 12/30/2012-END

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
				size_t uNrOfLocalCoefValues = 0;
				for(size_t h = 0; h < uNrOfLocalCoefs; h++)
					uNrOfLocalCoefValues += (size_t)pLocalCoefCounts[h];
				puStart[0] = (size_t)pLocalCoefOffsets[0];
				puCount[0] = uNrOfLocalCoefValues;

				// read the coefficients of the current basis
				#if	!WITH_SMART_PTR	// ADD-BY-LEETEN 12/30/2012
				TBuffer<TYPE_COEF_BIN> pLocalCoefBins;
				pLocalCoefBins.alloc(uNrOfLocalCoefValues);
				// ADD-BY-LEETEN 12/30/2012-BEGIN
				#else	// 	#if	!WITH_SMART_PTR
				boost::shared_array<TYPE_COEF_BIN> pLocalCoefBins(new TYPE_COEF_BIN[uNrOfLocalCoefValues]);
				#endif	// 	#if	!WITH_SMART_PTR
				// ADD-BY-LEETEN 12/30/2012-END

				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarCoefBin,
					puStart,
					puCount,
					(void*)&pLocalCoefBins[0]));

				#if	!WITH_SMART_PTR	// ADD-BY-LEETEN 12/30/2012
				TBuffer<TYPE_COEF_VALUE> pLocalCoefValues;
				pLocalCoefValues.alloc(uNrOfLocalCoefValues);
				// ADD-BY-LEETEN 12/30/2012-BEGIN
				#else	// #if	!WITH_SMART_PTR
				boost::shared_array<TYPE_COEF_VALUE> pLocalCoefValues(new TYPE_COEF_VALUE[uNrOfLocalCoefValues]);
				#endif	// #if	!WITH_SMART_PTR
				// ADD-BY-LEETEN 12/30/2012-END
				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarCoefValue,
					puStart,
					puCount,
					(void*)&pLocalCoefValues[0]));

				// scan through all basis
				for(size_t valuei = 0, lc = 0; lc < uNrOfLocalCoefs; lc++, i++)
				{
					vector<size_t> vuLocalCoefSub;
					_ConvertIndexToSub(lc, vuLocalCoefSub, vuLocalCoefLengths);

					// ADD-BY-LEETEN 12/26/2012-BEGIN
					size_t uCoefIndex = this->vuMapLocalToGlobal[i];
					bool bIsInCore = this->vbFlagsCoefInCore[uCoefIndex];

					if( !bIsInCore )
					{
						valuei += pLocalCoefCounts[lc];
						continue;
					}

					if( !vpcCoefPools[c] )
					{
						// MOD-BY-LEETEN 01/03/2013-FROM:						vpcCoefPools[c]= new CSepDWTPool<WT, unsigned short>;
						vpcCoefPools[c]= new CSepDWTPool<WT, BT>;
						// MOD-BY-LEETEN 01/03/2013-END
						vpcCoefPools[c]->_Set(
							// MOD-BY-LEETEN 01/03/2013-FROM:							UGetNrOfBins(),
							(BT)UGetNrOfBins(),
							// MOD-BY-LEETEN 01/03/2013-END
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
			LOG_VAR(uNrOfAllocatedPools);
			// ADD-BY-LEETEN 12/26/2012-END
			// #endif	// #if WITH_NETCDF 
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
						(void*)&pCoefBins[0]));
	
					ASSERT_NETCDF(nc_get_vara(
						iNcId,
						ncVarCoefValue,
						&uStart,
						&uCount,
						(void*)&pCoefValues[0]));
					for(size_t i = 0; i < uCount; i++)
						vpairCoefBinValues.push_back(pair<BT, WT>((BT)pCoefBins[i], (WT)pCoefValues[i]));
				}	
			}
		}
		// ADD-BY-LEETEN 01/05/2013-END

		//! Return the sum of all bins at the given position
		virtual	
		void
		_GetAllSums
		(
			const vector<size_t>& vuPos,
			// MOD-BY-LEETEN 01/03/2013-FROM:			vector<double>& vdSums,
			vector<ST>& vdSums,
			// MOD-BY-LEETEN 01/03/2013-END
			void *_Reserved = NULL
		)
		{
			vector<size_t> vuEmpty;	// this empty vector is used to speed up the access of elements
			// MOD-BY-LEETEN 01/05/2013-FROM:			vdSums.resize(UGetNrOfBins());
			vdSums.assign(UGetNrOfBins(), (ST)0);
			// MOD-BY-LEETEN 01/05/2013-END

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
			#if 0 			// MOD-BY-LEETEN 01/02/2013-FROM:
			uNrOfQueries++;
			#else // MOD-BY-LEETEN 01/02/2013-TO:
			this->uNrOfQueries++;
			#endif // MOD-BY-LEETEN 01/02/2013-END
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

				#if	0	// MOD-BY-LEETEN 01/05/2013-FROM:
				size_t uGlobalCoefIndex = UConvertSubToIndex(vuCoefSub, this->vuCoefLengths);
				if( this->vbFlagsCoefInCore[uGlobalCoefIndex] )
				{
					// MOD-BY-LEETEN 01/03/2013-FROM:					vector< pair<size_t, WT> > vpairCoefBinValues;
					vector< pair<BT, WT> > vpairCoefBinValues;
					// MOD-BY-LEETEN 01/03/2013-END
					this->vpcCoefPools[c]->_GetCoefSparse
					(
						vuLocalCoefSub,
						vpairCoefBinValues
					);

					// MOD-BY-LEETEN 01/03/2013-FROM:					for(vector< pair<size_t, double> >::iterator
					#if 0 // MOD-BY-LEETEN 01/04/2013-FROM:
					for(vector< pair<BT, WT> >::iterator
					      #else // MOD-BY-LEETEN 01/04/2013-TO:
					for(typename vector< pair<BT, WT> >::iterator
					      #endif // MOD-BY-LEETEN 01/04/2013-END
					// MOD-BY-LEETEN 01/03/2013-END
						ivpairCoefs = vpairCoefBinValues.begin();
						ivpairCoefs != vpairCoefBinValues.end();
						ivpairCoefs++ )
						vdSums[ivpairCoefs->first] += ivpairCoefs->second * dWavelet;
				}
				else
				{
					uNrOfIORequest++;	// ADD-BY-LEETEN 12/28/2012

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
				#else	// MOD-BY-LEETEN 01/05/2013-TO:
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
				#endif	// MOD-BY-LEETEN 01/05/2013-END
			}
			// ADD-BY-LEETEN 12/28/2012-BEGIN
			#if 0 			// MOD-BY-LEETEN 01/02/2013-FROM:
			uAccumNrOfIORequest += uNrOfIORequest;
			uMaxNrOfIORequest = max(uMaxNrOfIORequest, uNrOfIORequest);
			uMinNrOfIORequest = min(uMinNrOfIORequest, uNrOfIORequest);
			#else // MOD-BY-LEETEN 01/02/2013-TO:
			this->uAccumNrOfIORequest += uNrOfIORequest;
			this->uMaxNrOfIORequest = max(this->uMaxNrOfIORequest, uNrOfIORequest);
			#if	0	// MOD-BY-LEETEN 01/05/2012-FROM:
			this->uMinNrOfIORequest = min(this->uMinNrOfIORequest, uNrOfIORequest);
			#else	// ADD-BY-LEETEN 01/05/2012-TO:
			if( uNrOfIORequest )
			{
				if( !this->uMinNrOfIORequest )
					this->uMinNrOfIORequest = uNrOfIORequest;
				else
					this->uMinNrOfIORequest = min(this->uMinNrOfIORequest, uNrOfIORequest);
			}
			#endif	// ADD-BY-LEETEN 01/05/2012-END
			#endif // MOD-BY-LEETEN 01/02/2013-END
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
			#if	!WITH_PRE_ALLOCATED_QUEUES				// ADD-BY-LEETEN 01/05/2013
			vector< map< size_t, WT> > vmapLocalWavelet;
			vmapLocalWavelet.resize(uNrOfUpdatingCoefs);
			#endif	// #if	!WITH_PRE_ALLOCATED_QUEUES	// ADD-BY-LEETEN 01/05/2013

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

					#if	!WITH_PRE_ALLOCATED_QUEUES	// ADD-BY-LEETEN 01/05/2013
					typename map<size_t, WT>::iterator ivmapLocalWavelet = vmapLocalWavelet[c].find(uLocalCoef);
					if(vmapLocalWavelet[c].end() == ivmapLocalWavelet)
						vmapLocalWavelet[c].insert(pair<size_t, WT>(uLocalCoef, (WT)iSign * (WT)dWavelet));
					else
						ivmapLocalWavelet->second += (WT)iSign * (WT)dWavelet;
					// ADD-BY-LEETEN 01/05/2013-BEGIN
					#else	// #if	!WITH_PRE_ALLOCATED_QUEUES	
					this->vpairLocalWaveletQueues[c * uNrOfQueries + q] = pair<size_t, WT>(uLocalCoef, (WT)iSign * (WT)dWavelet);
					#endif	// #if	!WITH_PRE_ALLOCATED_QUEUES	
					// ADD-BY-LEETEN 01/05/2013-END
				}
			}

			#if	!WITH_PRE_ALLOCATED_QUEUES	// ADD-BY-LEETEN 01/05/2013
			for(size_t c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				for(typename map<size_t, WT>::iterator 
					ivmapLocalWavelet = vmapLocalWavelet[c].begin();
					ivmapLocalWavelet != vmapLocalWavelet[c].end();
					ivmapLocalWavelet++)
				{
					#if	0	// MOD-BY-LEETEN 01/05/2013-FROM:
					if( this->vpcCoefPools[c] )
					{
						vector< pair<BT, WT> > vpairCoefBinValues;
						this->vpcCoefPools[c]->_GetCoefSparse
						(
							ivmapLocalWavelet->first,
							vpairCoefBinValues
						);

						for(typename vector< pair<BT, WT> >::iterator
							ivpairCoefs = vpairCoefBinValues.begin();
							ivpairCoefs != vpairCoefBinValues.end();
							ivpairCoefs++ )
							vdSums[ivpairCoefs->first] += ivpairCoefs->second * ivmapLocalWavelet->second;
					}
					else
					{
						vector<size_t> vuLocalCoef;
						_ConvertIndexToSub(ivmapLocalWavelet->first, vuLocalCoef, this->vvuLocalLengths[c]);
						vector<size_t> vuGlobalCoef;
						vuGlobalCoef.resize(UGetNrOfDims());
						for(size_t d = 0; d < vuGlobalCoef.size(); d++)
							vuGlobalCoef[d] = this->vvuGlobalBase[c][d] + vuLocalCoef[d];
						size_t uGlobalCoefIndex = UConvertSubToIndex(vuGlobalCoef, this->vuCoefLengths);

						size_t uStart = (size_t)vuCoefOffsets[uGlobalCoefIndex];
						size_t uCount = (size_t)vusCoefCounts[uGlobalCoefIndex];

						if( uCount )
						{
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
						}	

						for(size_t i = 0; i < uCount; i++)
							vdSums[pCoefBins[i]] += pCoefValues[i] * ivmapLocalWavelet->second;
					}
					uNrOfCoefQueries++;
					#else	// MOD-BY-LEETEN 01/05/2013-TO:
					vector< pair<BT, WT> > vpairCoefBinValues;
					if(	this->dWaveletThreshold < fabs((double)ivmapLocalWavelet->second) )
					{
						this->_GetCoefSparse
						(
							c,
							ivmapLocalWavelet->first,
							vpairCoefBinValues
						);
						for(typename vector< pair<BT, WT> >::iterator
							ivpairCoefs = vpairCoefBinValues.begin();
							ivpairCoefs != vpairCoefBinValues.end();
							ivpairCoefs++ )
								vdSums[ivpairCoefs->first] += ivpairCoefs->second * ivmapLocalWavelet->second;
						uNrOfCoefQueries++;
					}
					#endif	// MOD-BY-LEETEN 01/05/2013-END
				}
			}
			// ADD-BY-LEETEN 01/05/2013-BEGIN
			#else	// #if	!WITH_PRE_ALLOCATED_QUEUES	
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
						vpairPrevLocalCoefQueues = vpairLocalWaveletQueues[i];
					}
				}
			}
			#endif	// #if	!WITH_PRE_ALLOCATED_QUEUES
			// ADD-BY-LEETEN 01/05/2013-END
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

		// ADD-BY-LEETEN 12/29/2012-BEGIN
		virtual
		void
		_DecodeBin
		(
			#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
			unsigned short usBin,
			valarray<DT> &vSAT,
			#else	// MOD-BY-LEETEN 01/03/2013-TO:
			const BT& usBin,
			valarray<ST> &vSAT,
			#endif	// MOD-BY-LEETEN 01/03/2013-END
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
				// if the coefficient of this wavelet is out of core, load them to the memory
				#if	!WITH_SMART_PTR				// ADD-BY-LEETEN 12/30/2012
				TBuffer<TYPE_COEF_BIN>	pCoefBins;		
				TBuffer<TYPE_COEF_VALUE> pCoefValues;	
				#endif	// #if	!WITH_SMART_PTR	// ADD-BY-LEETEN 12/30/2012
				if( !this->vpcCoefPools[w] )
				{
					size_t uLocalValueBase = this->vuGlobalValueBase[w];
					size_t uNrOfLocalValues = this->vuLocalValueCount[w];
					#if	!WITH_SMART_PTR	// ADD-BY-LEETEN 12/30/2012
					pCoefBins.alloc(uNrOfLocalValues);
					pCoefValues.alloc(uNrOfLocalValues);
					// ADD-BY-LEETEN 12/30/2012-BEGIN
					#else	// #if	!WITH_SMART_PTR
					boost::shared_array<TYPE_COEF_BIN> pCoefBins(new TYPE_COEF_BIN[uNrOfLocalValues]);
					boost::shared_array<TYPE_COEF_VALUE> pCoefValues(new TYPE_COEF_VALUE[uNrOfLocalValues]);
					#endif	// #if	!WITH_SMART_PTR
					// ADD-BY-LEETEN 12/30/2012-END
					ASSERT_NETCDF(nc_get_vara(
						iNcId,
						ncVarCoefBin,
						&uLocalValueBase,
						&uNrOfLocalValues,
						(void*)&pCoefBins[0]));
	
					ASSERT_NETCDF(nc_get_vara(
						iNcId,
						ncVarCoefValue,
						&uLocalValueBase,
						&uNrOfLocalValues,
						(void*)&pCoefValues[0]));

				}
				///////////////////////////////////////
				size_t uNrOfLocalCoefs = 1;
				_ConvertWaveletToLevels(w, vuGlobalCoefBase, vuLocalCoefLengths, uNrOfLocalCoefs);
				for(size_t lc = 0,		uLocalValueBase = 0, gc = 0;
					lc < uNrOfLocalCoefs; 
					lc++,				uLocalValueBase += (size_t)vusCoefCounts[gc], c++)
				{
					gc = vuMapLocalToGlobal[c];

					#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
					unsigned short usCount = vusCoefCounts[gc];
					unsigned short usNextOffset = vusCachedNextOffsets[gc];
					unsigned short usFetchedBin; // = vusCachedBins[gc];
					#else	// MOD-BY-LEETEN 01/03/2013-TO:
					BT usCount = vusCoefCounts[gc];
					BT usNextOffset = vusCachedNextOffsets[gc];
					BT usFetchedBin; // = vusCachedBins[gc];
					#endif	// MOD-BY-LEETEN 01/03/2013-END
					WT FetchedValue; // = vCachedValues[gc];

					// move the cache bin till it is great than or equal to the given bin
					if( !this->vbFlagsCoefInCore[gc] )
					{
						for(; usNextOffset < usCount; usNextOffset++)
						{
							// if( !usNextOffset || usBin < usFetchedBin)  
							usFetchedBin = pCoefBins[uLocalValueBase + (size_t)usNextOffset];
							FetchedValue = pCoefValues[uLocalValueBase + (size_t)usNextOffset];
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
			#if	0	// DEL-BY-LEETEN 01/05/2013-BEGIN
			vector<size_t> vuScanLineBase;
			vuScanLineBase.resize(UGetNrOfDims());

			vector<size_t> vuOtherCoefLengths;
			vuOtherCoefLengths.resize(UGetNrOfDims());
			#endif	// DEL-BY-LEETEN 01/05/2013-END
			      for(size_t uOffset = 1, d = 0, uCoefLength = 1; 
				d < UGetNrOfDims(); 
				uOffset *= uCoefLength, d++)
			{
			  uCoefLength = this->vuCoefLengths[d]; // ADD-BY-LEETEN 12/31/2012

				if( 1 == uCoefLength )
					continue;

				valarray<WT> vSrc;
				vSrc.resize(uCoefLength);

				valarray<WT> vDst;
				vDst.resize(uCoefLength);

				/*
				vector<size_t> vuScanLineIndices;
				vuScanLineIndices.resize(uCoefLength);
				*/

				size_t uNrOfLevels = vuDimLevels[d] - 1;
				size_t uNrOfScanLines = this->uNrOfCoefs / uCoefLength;

				for(size_t i = 0; i < uNrOfScanLines; i++)
				{
					size_t uScanlineBase = vvuSliceScanlineBase[d][i];
					vSrc = vSAT[slice(uScanlineBase, uCoefLength, uOffset)];
					_IDWT1D(vSrc, vDst, 2, uNrOfLevels - 1);
					vSAT[slice(uScanlineBase, uCoefLength, uOffset)] = vDst;
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
			#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
			const valarray<DT>& vCoefField,
			valarray<DT>& vDataField,
			#else	// MOD-BY-LEETEN 01/03/2013-TO:
			const valarray<ST>& vCoefField,
			valarray<ST>& vDataField,
			#endif	// MOD-BY-LEETEN 01/03/2013-END
			void* _Reserved = NULL
			)
		{
			// only keep the entropy field within the data range
			if( uDataSize != vDataField.size() )
				vDataField.resize(uDataSize);
			vector<size_t> vuSub;
			for(size_t d = 0; d < uDataSize; d++)
			{
				vector<size_t> vuSub;
				_ConvertIndexToSub(d, vuSub, vuDimLengths);
				vDataField[d] = vCoefField[UConvertSubToIndex(vuSub, vuCoefLengths)];
			}
		}

		virtual
		void
		_ClampBorder(
			// MOD-BY-LEETEN 01/03/2013-FROM:			valarray<DT>& vField,
			valarray<ST>& vField,
			// MOD-BY-LEETEN 01/03/2013-END
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
					// MOD-BY-LEETEN 01/03/2013-FROM:					vField[d] = (DT)0;
					vField[d] = (ST)0;
					// MOD-BY-LEETEN 01/03/2013-END
					continue;
				}
			}
		}
		// ADD-BY-LEETEN 01/02/2013-END

		#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
		CSATSepDWTDecoder()
			:bIsPrintingDecodeBinTiming(false)	// ADD-BY-LEETEN 12/29/2012-BEGIN
		#else	// MOD-BY-LEETEN 01/02/2013-TO:
		CSATSepDWTDecoder():
			// ADD-BY-LEETEN 01/05/2012-BEGIN
			uMinNrOfCoefQueries(0),
			uMaxNrOfCoefQueries(0),
			uAccumNrOfCoefQueries(0),
			uNrOfRangeQueries(0),
			// ADD-BY-LEETEN 01/05/2012-END
			#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
			CDecoderBase<DT>(),
			CSATSepDWTNetCDF(),
			CSepDWTHeader()
			#else	// MOD-BY-LEETEN 01/03/2013-TO:
			CDecoderBase<ST, BT>(),
			CSATSepDWTNetCDF(),
			// MOD-BY-LEETEN 01/04/2013-FROM:			CSepDWTHeader<ST, BT>()
			CSepDWTHeader()
			// MOD-BY-LEETEN 01/04/2013-END
			#endif	// MOD-BY-LEETEN 01/03/2013-END
		#endif	// MOD-BY-LEETEN 01/02/2013-END
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
