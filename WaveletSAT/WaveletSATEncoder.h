#pragma once

#include <vector>
using namespace std;
#include <math.h>

#include <string.h> 

#include "utils.h"		

#if		!defined(WITH_NETCDF)
#define WITH_NETCDF		1
#endif	// #if	!defined(WITH_NETCDF)
#if		!defined(WITH_NETCDF4)
#define WITH_NETCDF4	1
#endif	// #if	!defined(WITH_NETCDF4)

#include "SepDWTHeader.h"	
#include "SepDWTEncoderPool.h"
#include "EncoderBase.h"	
#include "SATSepDWTNetCDF.h"

#include "liblog.h"	

#if	!WITH_SMART_PTR	
#include "libbuf.h"

#else	// #if	!WITH_SMART_PTR	
#include <boost/shared_array.hpp>
#endif	// #if	!WITH_SMART_PTR	

#if	WITH_NETCDF
#include <netcdf.h>
#include "lognc.h"
#endif	// #if	WITH_NETCDF

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

	//! The WaveletSAT algorithm.
	/*!
	In order to apply wavelet transform for SATs or Integral Histograms, please follow the procedure below:
	1. Setup the data size by _SetDimLengths().
	2. Setup the #SATs by _AllocateBins().
	3. Call _Update() for each data point to incrementally compute wavelet.
	4. Call _Finalize() to finish the computation.

	To query the entry for a certain location:
	1. Call _GetAllSums() to get all SAT values for the given location.
	*/
	template<
		typename DT,	//!< Type of the data
		typename ST = typeSum,	//!< Type of the sum
		typename BT = typeBin,	//!< Type of the bin
		typename WT = typeWavelet	//!< Type of the wavelet coefficientsd
	>
	class CWaveletSATEncoder:
		virtual public CSATSepDWTNetCDF,
		virtual public CSepDWTHeader,
		virtual public CEncoderBase<DT, ST, BT>	
	{
protected:	
		size_t	uMinNrOfBufferedHeaders;
		size_t	uMaxMemoryUsage;
		#if		WITH_STREAMING		
		FILE *fpCoefBins;
		const char* szCoefBinTempFilename;
		FILE *fpCoefValues;
		const char* szCoefValueTempFilename;

		char szTempFilename[NC_MAX_NAME+1];

		size_t uNrOfWrittenCoefs;
		#endif	// #if		WITH_STREAMING		

	  //! Number of non-0 coefficients from all bin SATs.
	  size_t uNrOfNonZeroValues;

	  size_t uNrOfValuesInFullArray;

	  size_t uNrOfValuesOnSparseArray;

		//! The flag whether the method _Finalize() should multiple the result by the wavelet coefficients
		bool bIsFinalizedWithoutWavelet;

		//! #Coefs stored in full arrays
		size_t uNrOfCoefsInFullArray;	

		typedef CSepDWTEncoderPool<WT, BT> CSepDWTPool;
		vector< CSepDWTPool > vcCoefPools;

		virtual	
		void 
		_Update
		(
			const vector<size_t>& vuPos,
			const DT& value,
			void *_Reserved = NULL
		)
		{
			vector< pair<BT, ST> > vpBins;
			_MapValueToBins(vuPos, value, vpBins);
			if( vpBins.size() > 1 )
				for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
					this->vcCoefPools[c]._IncreaseCount(vpBins.size() - 1);

			for(typename vector< pair<BT, ST> >::iterator ivpBin  = vpBins.begin(); ivpBin != vpBins.end(); ivpBin++)
				_UpdateBin(vuPos, value, ivpBin->first, ivpBin->second);
		}

		#if	WITH_STREAMING
		virtual
		void
		_FlushBuffer
		(
			size_t c,
			void* _Reserved = NULL
		)
		{
			size_t uNrOfBufferedHeaders;
			size_t puStart[NC_MAX_DIMS];
			size_t puCount[NC_MAX_DIMS];
			TYPE_COEF_COUNT*	pCoefCounts;
			TYPE_COEF_OFFSET*	pCoefOffsets;
			size_t uNrOfBufferedCoefs;
			TYPE_COEF_VALUE*	pCoefValues;
			TYPE_COEF_BIN*		pCoefBins;

			this->vcCoefPools[c]._GetFileBuffer(
				puStart, 
				puCount, 
				uNrOfBufferedHeaders, 
				&pCoefCounts, 
				&pCoefOffsets, 
				uNrOfBufferedCoefs, 
				&pCoefValues, 
				&pCoefBins);

			for(size_t h = 0; h < uNrOfBufferedHeaders; h++)
				pCoefOffsets[h] += uNrOfWrittenCoefs;

			ASSERT_NETCDF(nc_put_vara(
								iNcId,
								ncVarCoefCount,
								puStart,
								puCount,
								(void*)&pCoefCounts[0] ));

			ASSERT_NETCDF(nc_put_vara(
								iNcId,
								ncVarCoefOffset,
								puStart,
								puCount,
								(void*)&pCoefOffsets[0] ));

			fwrite(&pCoefBins[0], sizeof(TYPE_COEF_BIN), uNrOfBufferedCoefs, fpCoefBins);
			fwrite(&pCoefValues[0], sizeof(TYPE_COEF_VALUE), uNrOfBufferedCoefs, fpCoefValues);
			uNrOfWrittenCoefs += uNrOfBufferedCoefs;
			size_t uMemoryUsage = WaveletSAT::UGetMemoryUsage();
			uMaxMemoryUsage = max(uMaxMemoryUsage, uMemoryUsage);
			this->vcCoefPools[c]._ResetFileBuffer();
		}
		#endif	// #if	WITH_STREAMING

		//! Update the specified bin.
		virtual 
		void 
		_UpdateBin
		(
			const vector<size_t>& vuPos, 
			const DT& value,
			const BT& uBin, 
			const ST& dWeight,
			void *_Reserved = NULL
		)
		{
			vector<size_t> vuLocalCoefSub;
			vuLocalCoefSub.resize(UGetNrOfDims());

			// now find the combination of the coefficients of all dimensions
			vdBinWeights[uBin] += dWeight;	

			vector<long> vlWavelets;
			vector<size_t> vuSubs;	// the subscripts of different level

			_GetForwardWavelet(vuPos, vuSubs, vlWavelets, true);

			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				long lWavelet = 1;
							
				for(size_t d = 0, uBase = 0;
					d < UGetNrOfDims(); 
					uBase += vuDimMaxLevels[d], d++, p++)
				{
					size_t uLevel = vuCoefDim2Level[p];
					lWavelet *= vlWavelets[uBase + uLevel];
					vuLocalCoefSub[d] = vuSubs[uBase + uLevel];
				}

				WT dWavelet = (WT)dWeight * (WT)lWavelet;
				this->vcCoefPools[c]._AddAt(uBin, vuLocalCoefSub, dWavelet);

				#if	WITH_STREAMING
				if( this->vcCoefPools[c].BIsReadyToFlush() ) 
					_FlushBuffer(c);
				#endif	//	#if	WITH_STREAMING
			}
		}
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0400,
			MIN_NR_OF_BUFFERED_HEADERS,	
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
		  CSATSepDWTNetCDF::_SetInteger(eName, lValue);
		  CSepDWTHeader::_SetInteger(eName, lValue);
		  // CEncoderBase<T, double>::_SetInteger(eName, lValue);
			
		  switch(eName) {
		  case MIN_NR_OF_BUFFERED_HEADERS:
			  uMinNrOfBufferedHeaders = (size_t)lValue;
			  if(uMinNrOfBufferedHeaders < 1) {
				  LOG("Warning: uMinNrOfBufferedHeaders is clamped to 1.");
				  uMinNrOfBufferedHeaders = 1;
			  }
			  break;
		  }
		}

		#if	WITH_STREAMING
		void
		_OpenFile
		(
			void *_Reserved = NULL
		)
		{
			fpCoefBins = fopen(szCoefBinTempFilename, "wb");
			fpCoefValues = fopen(szCoefValueTempFilename, "wb");
			#if WITH_NETCDF 
			#if !WITH_NETCDF4
			const char* szFileExt = "wnc";
			#else // #if !WITH_NETCDF4
			const char* szFileExt = "wnc4";
			#endif // #if !WITH_NETCDF4

			char szFilepath[NC_MAX_NAME+1];
			sprintf(szFilepath, "tmp.%s", szFileExt);
			strcpy(szTempFilename, szFilepath);
			// Create the file.
			#if !WITH_NETCDF4 
			ASSERT_NETCDF(nc_create(
    				szFilepath,
    				NC_CLOBBER,
    				&iNcId));
			#else	// #if !WITH_NETCDF4 
			ASSERT_NETCDF(nc_create(
    				szFilepath,
				NC_CLOBBER | NC_NETCDF4,
    				&iNcId));
			#endif	// #if !WITH_NETCDF4 

			// define the dimensions, including:
			// data_dim_0, data_dim_1, ... for the data
			// coef_dim_0, coef_dim_1, ... for the level
			// level_dim_0, level_dim_1, ... for the level
			size_t uDimLength;
			int iNcDimId;
			vncDims.clear();

			// define the dimension for #bins
			ASSERT_NETCDF(nc_def_dim(
						iNcId,
						szDimBin,
						(int)UGetNrOfBins(),
						&ncDimBin) );
			// define the dimension for #dims
			ASSERT_NETCDF(nc_def_dim(
						iNcId,
						szDimDim,
						(int)UGetNrOfDims(),
						&ncDimDim ) );

			for(size_t t = 0; t < NR_OF_DIM_TYPES; t++)
				for(int d = 0; d < UGetNrOfDims(); d++)
				{
					switch(t)
					{
					case DIM_TYPE_DATA:
						uDimLength = (size_t)this->vuDimLengths[d];
						break;
					case DIM_TYPE_LEVEL:
						uDimLength = (size_t)this->vuDimLevels[d];
						break;
					case DIM_TYPE_COEF:
						uDimLength = (size_t)((size_t)1<<(this->vuDimLevels[d]-1));
						break;
					}

					char szNcDimName[NC_MAX_NAME+1];
					sprintf(szNcDimName, "%s_DIM_%d", pszDimTypes[t], (unsigned int)d);

					ASSERT_NETCDF(nc_def_dim(
    								iNcId,
    								szNcDimName,
    								(int)uDimLength,
    								&iNcDimId));
					vncDims.push_back(iNcDimId);
				}

			//
			// now define the variable for the coef headers
			int piDimIds[NC_MAX_DIMS];
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				piDimIds[d] = vncDims[UGetNrOfDims() * DIM_TYPE_COEF + UGetNrOfDims() - 1 - d];

			ASSERT_NETCDF(nc_def_var(
					iNcId,
					szVarCoefCount,
					typeCoefCount,
					(int)UGetNrOfDims(),
					piDimIds,
					&ncVarCoefCount));

			ASSERT_NETCDF(nc_def_var(
					iNcId,
					szVarCoefOffset, 
					typeCoefOffset,
					(int)UGetNrOfDims(),
					piDimIds,
					&ncVarCoefOffset));

			#if WITH_NETCDF4 
			ASSERT_NETCDF(nc_def_var_deflate(
				   iNcId,
				   ncVarCoefCount, 
				   0, 
				   1, 
				   iDeflateLevel));
			ASSERT_NETCDF(nc_def_var_deflate(
				   iNcId,
				   ncVarCoefOffset, 
				   0, 
				   1, 
				   iDeflateLevel));
			#endif // #if WITH_NETCDF4

			// finish the definition mode
			ASSERT_NETCDF(nc_enddef(iNcId));
			#endif	// #if WITH_NETCDF 

			uNrOfWrittenCoefs = 0;
		}
		#endif	//	#if	WITH_STREAMING

		virtual 
		void
		_SaveFile
		(
			const char* szFilepathPrefix,
			void *_Reserved = NULL
		)
		{
			#if WITH_NETCDF 
			#if !WITH_NETCDF4
			const char* szFileExt = "wnc";
			#else // #if !WITH_NETCDF4
			const char* szFileExt = "wnc4";
			#endif // #if !WITH_NETCDF4

			char szFilepath[NC_MAX_NAME+1];
			sprintf(szFilepath, "%s.%s", szFilepathPrefix, szFileExt);
			#if	!WITH_STREAMING	
			// Create the file.
			#if !WITH_NETCDF4 
			ASSERT_NETCDF(nc_create(
    				szFilepath,
    				NC_CLOBBER,
    				&iNcId));
			#else	// #if !WITH_NETCDF4 
			ASSERT_NETCDF(nc_create(
    				szFilepath,
				NC_CLOBBER | NC_NETCDF4,
    				&iNcId));
			#endif	// #if !WITH_NETCDF4 

			// define the dimensions, including:
			// data_dim_0, data_dim_1, ... for the data
			// coef_dim_0, coef_dim_1, ... for the level
			// level_dim_0, level_dim_1, ... for the level
			size_t uDimLength;
			int iNcDimId;
			vncDims.clear();

			// define the dimension for coefficients, which is unlimited

			ASSERT_NETCDF(nc_def_dim(
						iNcId,
						szDimValue,
						uNrOfNonZeroValues,
						&ncDimValue) );

			// define the dimension for #bins
			ASSERT_NETCDF(nc_def_dim(
						iNcId,
						szDimBin,
						(int)UGetNrOfBins(),
						&ncDimBin) );
			// define the dimension for #bins
			ASSERT_NETCDF(nc_def_dim(
						iNcId,
						szDimDim,
						(int)UGetNrOfDims(),
						&ncDimDim ) );

			for(size_t t = 0; t < NR_OF_DIM_TYPES; t++)
				for(int d = 0; d < UGetNrOfDims(); d++)
				{
					switch(t)
					{
					case DIM_TYPE_DATA:
						uDimLength = (size_t)this->vuDimLengths[d];
						break;
					case DIM_TYPE_LEVEL:
						uDimLength = (size_t)this->vuDimLevels[d];
						break;
					case DIM_TYPE_COEF:
						uDimLength = (size_t)((size_t)1<<(this->vuDimLevels[d]-1));
						break;
					}

					char szNcDimName[NC_MAX_NAME+1];
					sprintf(szNcDimName, "%s_DIM_%d", pszDimTypes[t], (unsigned int)d);
					/*
					LOG_VAR(szNcDimName);
					LOG_VAR(uDimLength);
					*/

					ASSERT_NETCDF(nc_def_dim(
    								iNcId,
    								szNcDimName,
    								(int)uDimLength,
    								&iNcDimId));
					vncDims.push_back(iNcDimId);
				}

			//
			// now define the variable for the coef headers
			int piDimIds[NC_MAX_DIMS];
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				piDimIds[d] = vncDims[UGetNrOfDims() * DIM_TYPE_COEF + UGetNrOfDims() - 1 - d];

			ASSERT_NETCDF(nc_def_var(
					iNcId,
					szVarCoefCount,
					typeCoefCount,
					(int)UGetNrOfDims(),
					piDimIds,
					&ncVarCoefCount));

			ASSERT_NETCDF(nc_def_var(
					iNcId,
					szVarCoefOffset, 
					typeCoefOffset,
					(int)UGetNrOfDims(),
					piDimIds,
					&ncVarCoefOffset));

			// define the pool of the coefficients
			piDimIds[0] = this->ncDimValue;
			ASSERT_NETCDF(nc_def_var(
					iNcId,
					szVarCoefValue,
					typeCoefValue,
					1,
					piDimIds,
					&ncVarCoefValue));

			ASSERT_NETCDF(nc_def_var(
					iNcId,
					szVarCoefBin, 
					typeCoefBin,
					1,
					piDimIds,
					&ncVarCoefBin));

			#if WITH_NETCDF4 
			ASSERT_NETCDF(nc_def_var_deflate(
				   iNcId,
				   ncVarCoefCount, 
				   0, 
				   1, 
				   iDeflateLevel));
			ASSERT_NETCDF(nc_def_var_deflate(
				   iNcId,
				   ncVarCoefOffset, 
				   0, 
				   1, 
				   iDeflateLevel));
			ASSERT_NETCDF(nc_def_var_deflate(
				   iNcId,
				   ncVarCoefBin, 
				   0, 
				   1, 
				   iDeflateLevel));
			ASSERT_NETCDF(nc_def_var_deflate(
				   iNcId,
				   ncVarCoefValue, 
				   0, 
				   1, 
				   iDeflateLevel));
			#endif // #if WITH_NETCDF4

			// finish the definition mode
			ASSERT_NETCDF(nc_enddef(iNcId));

			// convert the basis id to its level subscript
			vector<size_t> vuLocalCoefLengths;
			vuLocalCoefLengths.resize(UGetNrOfDims());
			vector<size_t> vuGlobalCoefBase;
			vuGlobalCoefBase.resize(UGetNrOfDims());

			size_t puStart[NC_MAX_DIMS];
			size_t puCount[NC_MAX_DIMS];

			size_t uNrOfWrittenValues = 0;
			for(size_t c = 0; c < this->uNrOfUpdatingCoefs; c++)
			{
				vector<size_t> vuLevelSub;
				_ConvertIndexToSub(c, vuLevelSub, this->vuDimMaxLevels);

				size_t uNrOfLocalCoefs = 1;
				for(size_t d = 0; d < vuLevelSub.size(); d++)
				{
					// From this subscript, we can get the dim. length in this basis. 
					size_t uLevel = vuLevelSub[d];
					size_t uLen = (!uLevel)?1:(1<<(uLevel - 1));
					vuLocalCoefLengths[d] = uLen;
					uNrOfLocalCoefs *= uLen;

					// we can also decide it base in the n-dim. pool
					vuGlobalCoefBase[d] = (!uLevel)?0:(1<<(uLevel - 1));
				}

				vector< pair<BT, WT> > vpairLocalCoefBinValue;

				#if	!WITH_SMART_PTR
				TBuffer<TYPE_COEF_COUNT>	pLocalCoefCounts;
				TBuffer<TYPE_COEF_OFFSET>	pLocalCoefOffsets;
				pLocalCoefCounts.alloc(uNrOfLocalCoefs);
				pLocalCoefOffsets.alloc(uNrOfLocalCoefs);
				#else		// #if	!WITH_SMART_PTR
				boost::shared_array<TYPE_COEF_COUNT>	pLocalCoefCounts(new TYPE_COEF_COUNT[uNrOfLocalCoefs]);
				boost::shared_array<TYPE_COEF_OFFSET>	pLocalCoefOffsets(new TYPE_COEF_OFFSET[uNrOfLocalCoefs]);
				#endif		// #if	!WITH_SMART_PTR

				for(size_t bc = 0; bc < uNrOfLocalCoefs; bc++)
				{
					if( this->vcCoefPools[c].BIsEmpty(bc) )
						continue;

					vector<size_t> vuLocalCoefSub;
					_ConvertIndexToSub(bc, vuLocalCoefSub, vuLocalCoefLengths);

					vector< pair<BT, WT> > vpairCoefBinValue;

					this->vcCoefPools[c]._GetCoefSparse
					(
						vuLocalCoefSub,
						vpairCoefBinValue
					);

					pLocalCoefCounts[bc] = (TYPE_COEF_COUNT)vpairCoefBinValue.size();
					pLocalCoefOffsets[bc] = (TYPE_COEF_OFFSET)uNrOfWrittenValues + vpairLocalCoefBinValue.size();
					vpairLocalCoefBinValue.insert(vpairLocalCoefBinValue.end(), vpairCoefBinValue.begin(), vpairCoefBinValue.end());
				}

				// write the header
				for(size_t d = 0; d < vuLocalCoefLengths.size(); d++)
				  {
				    puStart[d] = vuGlobalCoefBase[UGetNrOfDims() - 1 - d];
				    puCount[d] = vuLocalCoefLengths[UGetNrOfDims() - 1 - d];
				  }

				ASSERT_NETCDF(nc_put_vara(
							      iNcId,
							      ncVarCoefCount,
							      puStart,
							      puCount,
							      (void*)&pLocalCoefCounts[0] ));

				ASSERT_NETCDF(nc_put_vara(
							      iNcId,
							      ncVarCoefOffset,
							      puStart,
							      puCount,
							      (void*)&pLocalCoefOffsets[0] ));
								  
				// write the coefficients
				puStart[0] = uNrOfWrittenValues;
				puCount[0] = vpairLocalCoefBinValue.size();

				#if	!WITH_SMART_PTR	
				TBuffer<TYPE_COEF_BIN> pLocalCoefBins;
				pLocalCoefBins.alloc(vpairLocalCoefBinValue.size());

				TBuffer<TYPE_COEF_VALUE> pLocalCoefValues;
				pLocalCoefValues.alloc(vpairLocalCoefBinValue.size());
				#else	// #if	!WITH_SMART_PTR	
				boost::shared_array<TYPE_COEF_BIN>		pLocalCoefBins(new TYPE_COEF_BIN[vpairLocalCoefBinValue.size()]);
				boost::shared_array<TYPE_COEF_VALUE>	pLocalCoefValues(new TYPE_COEF_VALUE[vpairLocalCoefBinValue.size()]);
				#endif	// #if	!WITH_SMART_PTR	

				for(size_t c = 0; c < vpairLocalCoefBinValue.size(); c++)
				  {
				    pLocalCoefBins[c] = (TYPE_COEF_BIN)vpairLocalCoefBinValue[c].first;
				    pLocalCoefValues[c] = (TYPE_COEF_VALUE)vpairLocalCoefBinValue[c].second;
				  }
				ASSERT_NETCDF(nc_put_vara(
							      iNcId,
							      ncVarCoefBin,
							      puStart,
							      puCount,
							      (void*)&pLocalCoefBins[0]));
								  
				ASSERT_NETCDF(nc_put_vara(
								 iNcId,
								 ncVarCoefValue,
								 puStart,
								 puCount,
								 (void*)&pLocalCoefValues[0]));
				uNrOfWrittenValues += vpairLocalCoefBinValue.size();
			}

			// close the file
			ASSERT_NETCDF(nc_close(iNcId));
			#else	// #if	!WITH_STREAMING	
			unlink(szFilepath);
			ASSERT_OR_LOG(-1 != rename(szTempFilename, szFilepath), perror(""));
			#endif	// #if	!WITH_STREAMING		
			#else	// #if WITH_NETCDF 
			#endif	// #if WITH_NETCDF 
			// write the data resolution
		}

		//! Finalize the computation of SAT
		virtual	
		void 
		_Finalize
		(
			void *_Reserved = NULL
		)
		{
			LOG_VAR(uMaxMemoryUsage);
			_ShowMemoryUsage(false); 

			for(size_t c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				vector<size_t> vuLevels;
				_ConvertIndexToSub(c, vuLevels, this->vuDimMaxLevels);

				this->vcCoefPools[c]._Finalize( (WT)1.0 );
				// this->vcCoefPools[c]._Weight( Wavelet / dWaveletDenomiator );
			}

			_ShowMemoryUsage(false);

			#if	!WITH_STREAMING		
			uNrOfNonZeroValues = 0;
			uNrOfValuesInFullArray = 0;
			uNrOfValuesOnSparseArray = 0;
			for(size_t c = 0; c < this->uNrOfUpdatingCoefs; c++)
			{
				size_t uF, uS;
				this->vcCoefPools[c]._GetArraySize(uF, uS, (WT)DGetThreshold());
				uNrOfValuesInFullArray += uF;
				uNrOfValuesOnSparseArray += uS;
			}
			uNrOfNonZeroValues = uNrOfValuesInFullArray + uNrOfValuesOnSparseArray;
			#else	//	#if	!WITH_STREAMING		
			fclose(fpCoefBins);
			fclose(fpCoefValues);

			ASSERT_NETCDF(nc_redef(iNcId));
			ASSERT_NETCDF(nc_def_dim(
						iNcId,
						szDimValue,
						uNrOfWrittenCoefs,
						&ncDimValue) );

			// define the pool of the coefficients
			int piDimIds[NC_MAX_DIMS];
			piDimIds[0] = this->ncDimValue;
			ASSERT_NETCDF(nc_def_var(
					iNcId,
					szVarCoefValue,
					typeCoefValue,
					1,
					piDimIds,
					&ncVarCoefValue));

			ASSERT_NETCDF(nc_def_var(
					iNcId,
					szVarCoefBin, 
					typeCoefBin,
					1,
					piDimIds,
					&ncVarCoefBin));

			#if WITH_NETCDF4 
			ASSERT_NETCDF(nc_def_var_deflate(
				   iNcId,
				   ncVarCoefBin, 
				   0, 
				   1, 
				   iDeflateLevel));
			ASSERT_NETCDF(nc_def_var_deflate(
				   iNcId,
				   ncVarCoefValue, 
				   0, 
				   1, 
				   iDeflateLevel));
			#endif	// #if WITH_NETCDF4 

			ASSERT_NETCDF(nc_enddef(iNcId));
			 
			size_t puStart[NC_MAX_DIMS];
			size_t puCount[NC_MAX_DIMS];
			LOG_VAR(uNrOfWrittenCoefs);

			const size_t uMaxNrOfCoefsToMove = 64000000;

			fpCoefBins = fopen(this->szCoefBinTempFilename, "rb");
			vector<TYPE_COEF_BIN> vBins;
			vBins.resize(uMaxNrOfCoefsToMove);
			for(size_t uNrOfMovedCoefs = 0; uNrOfMovedCoefs < uNrOfWrittenCoefs; uNrOfMovedCoefs += uMaxNrOfCoefsToMove) {
				puStart[0] = uNrOfMovedCoefs;
				puCount[0] = min(uNrOfWrittenCoefs - uNrOfMovedCoefs, uMaxNrOfCoefsToMove);
				fread(vBins.data(), sizeof(vBins[0]), puCount[0], fpCoefBins);
				ASSERT_NETCDF(nc_put_vara(
								iNcId,
								ncVarCoefBin,
								puStart,
								puCount,
								(void*)vBins.data()));
			}
			fclose(fpCoefBins);
			unlink(this->szCoefBinTempFilename);
			vBins.clear();

			fpCoefValues = fopen(this->szCoefValueTempFilename, "rb");
			vector<TYPE_COEF_VALUE> vValues;
			vValues.resize(uMaxNrOfCoefsToMove);
			for(size_t uNrOfMovedCoefs = 0; uNrOfMovedCoefs < uNrOfWrittenCoefs; uNrOfMovedCoefs += uMaxNrOfCoefsToMove) {
				puStart[0] = uNrOfMovedCoefs;
				puCount[0] = min(uNrOfWrittenCoefs - uNrOfMovedCoefs, uMaxNrOfCoefsToMove);
				fread(vValues.data(), sizeof(vValues[0]), puCount[0], fpCoefValues);
				ASSERT_NETCDF(nc_put_vara(
								iNcId,
								ncVarCoefValue,
								puStart,
								puCount,
								(void*)vValues.data()));
			}
			fclose(fpCoefValues);
			unlink(this->szCoefValueTempFilename);
			vValues.clear();

			ASSERT_NETCDF(nc_close(iNcId));
			#endif	//	#if	!WITH_STREAMING		
		}

		//! Compute and display statistics for the computed wavelet coefficients.
		virtual	
		void
		_ShowStatistics
		(
			void *_Reserved = NULL
		)
		{
			LOG_VAR(uNrOfValuesInFullArray);
			LOG_VAR(uNrOfValuesOnSparseArray);
			LOG_VAR(uNrOfNonZeroValues);

			size_t uNrOfDataItems = this->uDataSize;

			double dCR = (double)(uNrOfDataItems * UGetNrOfBins()) / (double)uNrOfNonZeroValues;
			LOG_VAR(dCR);

			double dOverhead = (double)uNrOfNonZeroValues / (double)uNrOfDataItems;
			LOG_VAR(dOverhead);
		}

		//! Allocate the space to store coefficients for all bins. 
		/*! 
		*/
		virtual	
		void 
		_Allocate
		(
			void *_Reserved = NULL
		)
		{
			// compute the #coef in the full array per bin SAT
			uNrOfCoefsInFullArray = min(
				(size_t)floor( (double)uSizeOfFullArrays/(double)(UGetNrOfBins() * sizeof(WT))), 
				uNrOfCoefs);

			////////////////////////////////////////
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
				double dCurrentAspectRatio = -1.0; 
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

			///////////////////////////////
			vector<size_t> vuPoolDimLengths;
			vuPoolDimLengths.resize(this->UGetNrOfDims());
			this->vcCoefPools.resize(this->uNrOfUpdatingCoefs);
			for(size_t c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				vector<size_t> vuLocalCoefSub;
				_ConvertIndexToSub(c, vuLocalCoefSub, this->vuDimLevels);

				bool bIsSparse = false;
				// decide the pool size
				for(size_t d = 0; d < this->UGetNrOfDims(); d++)
				{
					vuPoolDimLengths[d] = (!vuLocalCoefSub[d])?1:((size_t)1 << (vuLocalCoefSub[d] - 1));
					// decide whether the array is sparse.
					if( vuLocalCoefSub[d] > vuOptimalDimLevel[d] )
						bIsSparse = true;
				}

				if( !uNrOfCoefsInFullArray )
				  bIsSparse = true;

				// decide the wavelet weight
				double dWavelet = +1.0;

				vector<size_t> vuWaveletLengths;
				vuWaveletLengths.resize(UGetNrOfDims());

				vector<size_t> vuGlobalCoefBase;

				for(size_t d = 0; d < vuLocalCoefSub.size(); d++)
				{
					size_t uLevel = vuLocalCoefSub[d];
					if( uLevel >= 1 )
						dWavelet *= (double)(1 << (uLevel - 1));

					vuWaveletLengths[d] = (size_t)1<<(( !uLevel )?(vuDimLevels[d] - 1):(vuDimLevels[d] - uLevel));

					#if		WITH_STREAMING		
					vuGlobalCoefBase.push_back( (!uLevel)?0:(1<<(uLevel - 1)) );
					#endif	//	#if	WITH_STREAMING		
				}

				dWavelet = sqrt(dWavelet);
				double dWeight = dWavelet/dWaveletDenomiator;
				this->vcCoefPools[c]._SetWaveletWeight(dWeight);
				this->vcCoefPools[c]._SetDataDimLengths(vuDimLengths);
				this->vcCoefPools[c]._SetWaveletLengths(vuWaveletLengths);
				#if	WITH_DYNAMIC_ARRAY_ALLOCATION		
				#if		!WITH_STREAMING		

				#else	// #if		!WITH_STREAMING		
				this->vcCoefPools[c]._SetBuffer
				(
					uMinNrOfBufferedHeaders,	
					vuGlobalCoefBase
				);
				#endif	// #if		!WITH_STREAMING		

				bIsSparse = true;
				#endif	// #if	WITH_DYNAMIC_ARRAY_ALLOCATION		

				this->vcCoefPools[c]._Set(
					(BT)this->UGetNrOfBins(),
					vuPoolDimLengths,
					this->vuMaxCounts[c],
					bIsSparse);
			}

			_ShowMemoryUsage(false); 

			#if		WITH_STREAMING		
			_OpenFile();
			#endif	// #if		WITH_STREAMING		
		}
		
		CWaveletSATEncoder():
			uMaxMemoryUsage(0),	
			uMinNrOfBufferedHeaders(8),
			#if		WITH_STREAMING		
			szCoefBinTempFilename("coef.bin.tmp"),
			szCoefValueTempFilename("coef.value.tmp"),
			uNrOfWrittenCoefs(0),	
			#endif	// #if		WITH_STREAMING		
			bIsFinalizedWithoutWavelet(false)
		{
		}
	};
}

/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
