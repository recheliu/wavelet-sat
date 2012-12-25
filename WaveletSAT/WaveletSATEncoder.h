#pragma once

#define WITH_VECTORS_FOR_COUNTED_COEFS	0

// ADD-BY-LEETEN 10/10/2012-BEGIN
#define	WITH_BOUNDARY_AWARE_DWT		0
// ADD-BY-LEETEN 10/10/2012-END

// ADD-BY-LEETEN 11/11/2012-BEGIN
#define WITH_COEF_POOL			1
// ADD-BY-LEETEN 11/11/2012-END

#include <map>	

#include <vector>
using namespace std;
#include <math.h>

#include "utils.h"		// ADD-BY-LEETEN 10/29/2012
#include "SepDWTHeader.h"	
#include "SepDWTData.h"		// ADD-BY-LEETEN 10/29/2012
#include "SepDWTPool.h"		// ADD-BY-LEETEN 11/11/2012
#include "EncoderBase.h"	
#include "SATSepDWTNetCDF.h"	// ADD-BY-LEETEN 12/16/2012

#include "liblog.h"	

// ADD-BY-LEETEN 12/12/2012-BEGIN
#include "libbuf.h"
#if	WITH_NETCDF
#include <netcdf.h>
#include "lognc.h"
#endif	// #if	WITH_NETCDF
// ADD-BY-LEETEN 12/12/2012-END

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
	template<typename T>
	class CWaveletSATEncoder:
		public CSATSepDWTNetCDF,	// ADD-BY-LEETEN 12/16/2012
		public CSepDWTHeader,
		public CEncoderBase<T, double>	
	{
protected:	
	  // ADD-BY-LEETEN 12/15/2012-BEGIN
	  //! Number of non-0 coefficients from all bin SATs.
	  size_t uNrOfNonZeroCoefs;

	  size_t uCountInFullArray;

	  size_t uCountInSparseArray;
	  // ADD-BY-LEETEN 12/15/2012-END

		//! The flag whether the method _Finalize() should multiple the result by the wavelet coefficients
		bool bIsFinalizedWithoutWavelet;

		//! #Coefs stored in full arrays
		size_t uNrOfCoefsInFullArray;	

		// ADD-BY-LEETEN 10/29/2012-BEGIN
		#if	!WITH_COEF_POOL	// ADD-BY-LEETEN 11/11/2012
		vector< CSepDWTData<double> > vcBinCoefs;
		// ADD-BY-LEETEN 11/11/2012-BEGIN
		#else	// #if	!WITH_COEF_POOL	
		vector< CSepDWTPool<double, unsigned short> > vcCoefPools;
		#endif	// #if	!WITH_COEF_POOL	
		// ADD-BY-LEETEN 11/11/2012-END
		// ADD-BY-LEETEN 10/29/2012-END

		//! Update the specified bin.
		void 
		_UpdateBin
		(
			const vector<size_t>& vuPos, 
			const T& value,
			size_t uBin, 
			const double& dWeight,
			void *_Reserved = NULL
		)
		{
			#if	!WITH_COEF_POOL	// ADD-BY-LEETEN 11/11/2012
			// ADD-BY-LEETEN 10/29/2012-BEGIN
			#if		!WITH_1D_DIVISION
			vector<size_t> vuCoefPos;
			vuCoefPos.resize(UGetNrOfDims());
			#endif	// #if	!WITH_1D_DIVISION
			// ADD-BY-LEETEN 10/29/2012-END

			vdBinWeights[uBin] += dWeight;	// ADD-BY-LEETEN 10/10/2012

			#if	!WITH_PRECOMPUTED_WAVELET_SUMS	// ADD-BY-LEETEN 10/21/2012
			vector<long> vlWavelets;

			vlWavelets.resize(UGetNrOfDims() * this->vuDimLevels[0]);	// ADD-BY-LEETEN 10/19/2012

			// for each dimension, fetch the l[d] indices;
			for(size_t p = 0, d = 0; d < UGetNrOfDims(); d++)
			{
				size_t uDimMaxLevel = vuDimMaxLevels[d];
				size_t uPos = vuPos[d];
				size_t uMaxWin = 1 << vuDimLevels[d];
				for(size_t 	
					l = 0, w = uMaxWin;
					l < uDimMaxLevel; 
					l++, p++, w >>= 1)
				{
					long lWavelet;
					// Given the currenet level l and subscript uPos, compute the sum of the portion in wavelet after uPos
					size_t uPosInWavelet = uPos % w;
					if( 0 == l )
						lWavelet = (long)w / 2 - (long)uPosInWavelet;
					else
					{
						if( uPosInWavelet < w / 2)
							lWavelet = (long)uPosInWavelet;
						else
							lWavelet = (long)(w - uPosInWavelet);
						lWavelet *= -1;		
					}
					vlWavelets[p] = lWavelet;
				}
			}

			// now find the combination of the coefficients of all dimensions
			vector<size_t> vuPosLevelProduct;
			vuPosLevelProduct.resize(UGetNrOfDims());
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				vuPosLevelProduct[d] = vuPos[d] * vuDimLevels[d];

			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				long lWavelet = 1;
							
				#if	WITH_1D_DIVISION		
				size_t uCoefId = 0;
				#endif	// #if	WITH_1D_DIVISION	

				#if	!WITH_COEF_DIM_2_WAVELET	// ADD-BY-LEETEN 10/21/2012
				for(size_t d = 0, uBase = 0;
					d < UGetNrOfDims(); 
					uBase += vuDimMaxLevels[d], d++, p++)
				#else	// #if	!WITH_COEF_DIM_2_WAVELET	
				for(size_t d = 0; d < UGetNrOfDims(); d++, p++)
				#endif	// #if	!WITH_COEF_DIM_2_WAVELET	
				{
					size_t uLevel = vuCoefDim2Level[p];
					
					/*
					// skip the level if it is larger than the threshold 
					if( uLevel >= vuDimMaxLevels[d] )
						continue;
					*/
					
					#if	!WITH_COEF_DIM_2_WAVELET	// ADD-BY-LEETEN 10/21/2012
					// ADD-BY-LEETEN 10/18/2012-BEGIN
					if( 0 == d )
						lWavelet = vlWavelets[uBase + uLevel];
					else
					// ADD-BY-LEETEN 10/18/2012-END
					lWavelet *= vlWavelets[uBase + uLevel];
					// ADD-BY-LEETEN 10/21/2012-BEGIN
					#else	// #if	!WITH_COEF_DIM_2_WAVELET	
					size_t uWavelet = this->vuCoefDim2Wavelet[p];
					if( 0 == d )
						lWavelet = vlWavelets[uWavelet];
					else
						lWavelet *= vlWavelets[uWavelet];
					#endif	// #if	!WITH_COEF_DIM_2_WAVELET	
					// ADD-BY-LEETEN 10/21/2012-END

					#if	WITH_1D_DIVISION
					uCoefId += vvuSubLevel2Coef[d][vuPosLevelProduct[d] + uLevel];

					// ADD-BY-LEETEN 10/29/2012-BEGIN
					#else	// #if	WITH_1D_DIVISION
					vuCoefPos[d] = vvuSubLevel2Coef[d][vuPosLevelProduct[d] + uLevel];
					#endif	// #if	WITH_1D_DIVISION
					// ADD-BY-LEETEN 10/29/2012-END

					// ADD-BY-LEETEN 10/18/2012-BEGIN
					if( !lWavelet )
					{
						p += UGetNrOfDims() - d;
						break;
					}
					// ADD-BY-LEETEN 10/18/2012-END
				}
				if( !lWavelet )
					continue;
				double dWavelet = dWeight * (double)lWavelet;
			// ADD-BY-LEETEN 10/21/2012-BEGIN
			#else	// #if	!WITH_PRECOMPUTED_WAVELET_SUMS	
			vector< size_t > vuWaveletSums;
			vuWaveletSums.resize(UGetNrOfDims() * this->vuDimLevels[0]);

			// for each dimension, fetch the l[d] indices;
			for(size_t p = 0, uBase = 1, d = 0; d < UGetNrOfDims(); uBase *= (this->vuDimLengths[d] + 1), d++)
				for(size_t l = 0, w = 1 << vuDimLevels[d]; l < vuDimMaxLevels[d]; l++, w >>= 1, p++)
				{
					size_t uWaveletSum;
					// Given the currenet level l and subscript uPos, compute the sum of the portion in wavelet after uPos
					size_t uPosInWavelet = vuPos[d] % w;
					if( 0 == l )
						uWaveletSum = (w/2 - uPosInWavelet);
					else
						if( uPosInWavelet < w/2 )
							uWaveletSum = uPosInWavelet;
						else
							uWaveletSum = w - uPosInWavelet;

					vuWaveletSums[p] = uWaveletSum * uBase;
				}

			// now find the combination of the coefficients of all dimensions
			vector<size_t> vuPosLevelProduct;
			vuPosLevelProduct.resize(UGetNrOfDims());
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				vuPosLevelProduct[d] = vuPos[d] * vuDimLevels[d];

			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				#if	WITH_1D_DIVISION		
				size_t uCoefId = 0;
				#endif	// #if	WITH_1D_DIVISION	
				size_t uWaveletSumId = 0;
				#if	!WITH_COEF_DIM_2_WAVELET	// ADD-BY-LEETEN 10/21/2012-BEGIN
				for(size_t d = 0, uBase = 0; d < UGetNrOfDims(); uBase += vuDimMaxLevels[d], d++, p++)
				// ADD-BY-LEETEN 10/21/2012-BEGIN
				#else	// #if	!WITH_COEF_DIM_2_WAVELET	
				for(size_t d = 0; d < UGetNrOfDims(); d++, p++)
				#endif	// #if	!WITH_COEF_DIM_2_WAVELET	
				// ADD-BY-LEETEN 10/21/2012-END
				{
					size_t uLevel = vuCoefDim2Level[p];
					#if	WITH_1D_DIVISION	
					uCoefId += vvuSubLevel2Coef[d][vuPosLevelProduct[d] + uLevel];
					
					// ADD-BY-LEETEN 10/29/2012-BEGIN
					#else	// #if	WITH_1D_DIVISION
					vuCoefPos[d] = vvuSubLevel2Coef[d][vuPosLevelProduct[d] + uLevel]
					#endif	// #if	WITH_1D_DIVISION
					// ADD-BY-LEETEN 10/29/2012-END

					#if	!WITH_COEF_DIM_2_WAVELET	// ADD-BY-LEETEN 10/21/2012
					uWaveletSumId += vuWaveletSums[uBase + uLevel];
					// ADD-BY-LEETEN 10/21/2012-BEGIN
					#else	// #if	!WITH_COEF_DIM_2_WAVELET
					uWaveletSumId += vuWaveletSums[vuCoefDim2Wavelet[p]];
					#endif	// #if	!WITH_COEF_DIM_2_WAVELET
					// ADD-BY-LEETEN 10/21/2012-END
				}
				double dWavelet = this->vdWaveletSums[uWaveletSumId];
				dWavelet *= dWeight * this->vdWaveletSigns[c];
			#endif	// #if	!WITH_PRECOMPUTED_WAVELET_SUMS	
			// ADD-BY-LEETEN 10/21/2012-END

				// ADD-BY-LEETEN 10/31/2012-BEGIN
				#if	WITH_1D_DIVISION
				this->vcBinCoefs[uBin]._AddAtIndex(uCoefId, dWavelet);
				// ADD-BY-LEETEN 10/31/2012-END
				#else	// #if	WITH_1D_DIVISION
				this->vcBinCoefs[uBin]._AddAtPos(vuCoefPos, dWavelet);
				#endif	// #if	WITH_1D_DIVISION
				// ADD-BY-LEETEN 10/29/2012-END
			}
			// ADD-BY-LEETEN 11/11/2012-BEGIN
			#else	// #if	!WITH_COEF_POOL	
			vector<size_t> vuPoolSubs;
			vuPoolSubs.resize(UGetNrOfDims());

			// now find the combination of the coefficients of all dimensions
			vdBinWeights[uBin] += dWeight;	// ADD-BY-LEETEN 10/10/2012

			vector<long> vlWavelets;
			vlWavelets.resize(uNrOfWaveletsFromAllDims);

			vector<size_t> vuSubs;	// the subscripts of different level
			vuSubs.resize(uNrOfWaveletsFromAllDims);

			// for each dimension, fetch the l[d] indices;
			for(size_t p = 0, d = 0; d < UGetNrOfDims(); d++)
			{
				size_t uDimMaxLevel = vuDimMaxLevels[d];
				size_t uPos = vuPos[d];
				size_t uMaxWin = 1 << vuDimLevels[d];
				for(size_t 	
					l = 0, w = uMaxWin;
					l < uDimMaxLevel; 
					l++, p++, w >>= 1)
				{
					long lWavelet;
					// Given the currenet level l and subscript uPos, compute the sum of the portion in wavelet after uPos
					size_t uPosInWavelet = uPos % w;
					if( 0 == l )
						lWavelet = (long)w / 2 - (long)uPosInWavelet;
					else
					{
						if( uPosInWavelet < w / 2)
							lWavelet = (long)uPosInWavelet;
						else
							lWavelet = (long)(w - uPosInWavelet);
						lWavelet *= -1;		
					}
					vlWavelets[p] = lWavelet;
					vuSubs[p] = uPos / w;
				}
			}

			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				long lWavelet = 1;
							
				for(size_t d = 0, uBase = 0;
					d < UGetNrOfDims(); 
					uBase += vuDimMaxLevels[d], d++, p++)
				{
					size_t uLevel = vuCoefDim2Level[p];
					lWavelet *= vlWavelets[uBase + uLevel];
					vuPoolSubs[d] = vuSubs[uBase + uLevel];
				}

				double dWavelet = dWeight * (double)lWavelet;
				
				this->vcCoefPools[c]._AddAt(uBin, vuPoolSubs, dWavelet);
			}
			#endif	// #if	!WITH_COEF_POOL	
			// ADD-BY-LEETEN 11/11/2012-END
		}
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0400,
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
		  // CEncoderBase<T, double>::_SetLong(eName, lValue);
		}

		virtual	
		void
		_SetBoolean(
			int eName,
			bool bValue,
			void* _Reserved = NULL
		)
		{
		}

		// ADD-BY-LEETEN 12/12/2012-BEGIN
		virtual 
		void
		_SaveFile
		(
			const char* szFilepathPrefix,
			void *_Reserved = NULL
		)
		{
			#if WITH_NETCDF 
			// ADD-BY-LEETEN 12/15/2012-BEGIN
			#if !WITH_NETCDF4
			const char* szFileExt = "wnc";
			#else // #if !WITH_NETCDF4
			const char* szFileExt = "wnc4";
			#endif // #if !WITH_NETCDF4

			char szFilepath[NC_MAX_NAME+1];
			sprintf(szFilepath, "%s.%s", szFilepathPrefix, szFileExt);
			// ADD-BY-LEETEN 12/15/2012-END

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
			char* szNcDimName;
			int iNcDimId;
			vncDims.clear();

			// define the dimension for coefficients, which is unlimited

			ASSERT_NETCDF(nc_def_dim(
						iNcId,
						szDimValue,
						uNrOfNonZeroCoefs,
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
						uDimLength = (int)this->vuDimLengths[d];
						break;
					case DIM_TYPE_LEVEL:
						uDimLength = (int)this->vuDimLevels[d];
						break;
					case DIM_TYPE_COEF:
						uDimLength = (int)1<<(this->vuDimLevels[d]-1);
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

			#if	0	// DEL-BY-LEETEN 12/25/2012-BEGIN
			#if !WITH_NETCDF4
			typeHeaderOffset= NC_INT;
			typeHeaderCount	= NC_INT;
			typeCoefBin		= NC_INT;	
			#else // #if !WITH_NETCDF4
			typeHeaderOffset= NC_UINT64;
			typeHeaderCount	= NC_UINT;
			typeCoefBin		= NC_UINT;	
			#endif // #if !WITH_NETCDF4
			typeCoefValue	= NC_DOUBLE;
			#endif	// DEL-BY-LEETEN 12/25/2012-END
			ASSERT_NETCDF(nc_def_var(
					iNcId,
					szVarHeaderCount,
					typeHeaderCount,
					(int)UGetNrOfDims(),
					piDimIds,
					&ncVarHeaderCount));

			ASSERT_NETCDF(nc_def_var(
					iNcId,
					szVarHeaderOffset, 
					typeHeaderOffset,
					(int)UGetNrOfDims(),
					piDimIds,
					&ncVarHeaderOffset));

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

			// ADD-BY-LEETEN 12/16/2012-BEGIN
			#if WITH_NETCDF4 
			ASSERT_NETCDF(nc_def_var_deflate(
				   iNcId,
				   ncVarHeaderCount, 
				   0, 
				   1, 
				   iDeflateLevel));
			ASSERT_NETCDF(nc_def_var_deflate(
				   iNcId,
				   ncVarHeaderOffset, 
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
			// ADD-BY-LEETEN 12/16/2012-END

			// finish the definition mode
			ASSERT_NETCDF(nc_enddef(iNcId));

			// convert the basis id to its level subscript
			vector<size_t> vuBasisCoefLengths;
			vuBasisCoefLengths.resize(UGetNrOfDims());
			vector<size_t> vuBasisCoefBase;
			vuBasisCoefBase.resize(UGetNrOfDims());

			size_t puStart[NC_MAX_DIMS];
			size_t puCount[NC_MAX_DIMS];

			size_t uCoefBase = 0;
			for(size_t c = 0; c < this->uNrOfUpdatingCoefs; c++)
			{
				vector<size_t> vuLevelSub;
				_ConvertIndexToSub(c, vuLevelSub, this->vuDimMaxLevels);

				size_t uNrOfBasisCoefs = 1;
				for(size_t d = 0; d < vuLevelSub.size(); d++)
				{
					// From this subscript, we can get the dim. length in this basis. 
					size_t uLevel = vuLevelSub[d];
					size_t uLen = (!uLevel)?1:(1<<(uLevel - 1));
					vuBasisCoefLengths[d] = uLen;
					uNrOfBasisCoefs *= uLen;

					// we can also decide it base in the n-dim. pool
					vuBasisCoefBase[d] = (!uLevel)?0:(1<<(uLevel - 1));
				}

				vector< pair<size_t, double> > vpairCoefsInBasis;

				#if	0	// DEL-BY-LEETEN 12/25/2012-BEGIN
				#if !WITH_NETCDF4
				typedef int	TYPE_HEADER_COUNT;
				typedef int	TYPE_HEADER_OFFSET;
				typedef int	TYPE_COEF_BIN;
				#else // #if !WITH_NETCDF4
				typedef unsigned int 		TYPE_HEADER_COUNT;
				typedef unsigned long long	TYPE_HEADER_OFFSET;
				typedef unsigned int		TYPE_COEF_BIN;
				#endif // #if !WITH_NETCDF4
				typedef double TYPE_COEF_VALUE;
				#endif	// DEL-BY-LEETEN 12/25/2012-END
				TBuffer<TYPE_HEADER_COUNT>	pHeaderCount;
				TBuffer<TYPE_HEADER_OFFSET>	pHeaderOffset;
				
				pHeaderCount.alloc(uNrOfBasisCoefs);
				pHeaderOffset.alloc(uNrOfBasisCoefs);

				for(size_t bc = 0; bc < uNrOfBasisCoefs; bc++)
				{
					vector<size_t> vuBasisCoefSub;
					_ConvertIndexToSub(bc, vuBasisCoefSub, vuBasisCoefLengths);

					vector< pair<size_t, double> > vpairCoefs;
					this->vcCoefPools[c]._GetCoefSparse
					(
						vuBasisCoefSub,
						vpairCoefs
					);

					pHeaderCount[bc] = (TYPE_HEADER_COUNT)vpairCoefs.size();
					// MOD-BY-LEETEN 12/25/2012-FROM:	pHeaderOffset[bc] = (TYPE_HEADER_OFFSET)uCoefBase;
					pHeaderOffset[bc] = (TYPE_HEADER_OFFSET)uCoefBase + vpairCoefsInBasis.size();
					// MOD-BY-LEETEN 12/25/2012-END
					vpairCoefsInBasis.insert(vpairCoefsInBasis.end(), vpairCoefs.begin(), vpairCoefs.end());
				}
				// ADD-BY-LEETEN 12/15/2012-BEGIN
				// write the header
				for(size_t d = 0; d < vuBasisCoefLengths.size(); d++)
				  {
				    puStart[d] = vuBasisCoefBase[UGetNrOfDims() - 1 - d];
				    puCount[d] = vuBasisCoefLengths[UGetNrOfDims() - 1 - d];
				  }

				ASSERT_NETCDF(nc_put_vara(
							      iNcId,
							      ncVarHeaderCount,
							      puStart,
							      puCount,
							      (void*)&pHeaderCount[0] ));

				ASSERT_NETCDF(nc_put_vara(
							      iNcId,
							      ncVarHeaderOffset,
							      puStart,
							      puCount,
							      (void*)&pHeaderOffset[0] ));
								  
				// write the coefficients
				puStart[0] = uCoefBase;
				puCount[0] = vpairCoefsInBasis.size();

				TBuffer<TYPE_COEF_BIN> pCoefBins;
				pCoefBins.alloc(vpairCoefsInBasis.size());

				TBuffer<TYPE_COEF_VALUE> pCoefValues;
				pCoefValues.alloc(vpairCoefsInBasis.size());

				for(size_t c = 0; c < vpairCoefsInBasis.size(); c++)
				  {
				    pCoefBins[c] = (TYPE_COEF_BIN)vpairCoefsInBasis[c].first;
				    pCoefValues[c] = (TYPE_COEF_VALUE)vpairCoefsInBasis[c].second;
				  }

				ASSERT_NETCDF(nc_put_vara(
							      iNcId,
							      ncVarCoefBin,
							      puStart,
							      puCount,
							      (void*)&pCoefBins[0]));
								  
				ASSERT_NETCDF(nc_put_vara(
								 iNcId,
								 ncVarCoefValue,
								 puStart,
								 puCount,
								 (void*)&pCoefValues[0]));

				uCoefBase += vpairCoefsInBasis.size();
				// ADD-BY-LEETEN 12/15/2012-END
			}

			// close the file
			ASSERT_NETCDF(nc_close(iNcId));
			#else	// #if WITH_NETCDF 
			#endif	// #if WITH_NETCDF 
			// write the data resolution

		}
		// ADD-BY-LEETEN 12/12/2012-END

		//! Finalize the computation of SAT
		virtual	
		void 
		_Finalize
		(
			void *_Reserved = NULL
		)
		{
			
			_ShowMemoryUsage(false); // ADD-BY-LEETEN 11/14/2012

			#if	!WITH_COEF_POOL	// ADD-BY-LEETEN 11/11/2012
			vector<size_t> vuSub;	// ADD-BY-LEETEN 10/06/2012

			// ADD-BY-LEETEN 10/10/2012-BEGIN
			#if	WITH_BOUNDARY_AWARE_DWT	
			// consider the out-of-bound case
			vector<long> vlWavelets;
			// for each dimension, fetch the l[d] indices;
			for(size_t p = 0, d = 0; d < UGetNrOfDims(); d++)
			{
				size_t uDimMaxLevel = vuDimMaxLevels[d];
				size_t uPos = this->vuDataDimLengths[d] - 1;
				size_t uMaxWin = 1 << vuDimLevels[d];
				for(size_t 	
					l = 0, w = uMaxWin;
					l < uDimMaxLevel; 
					l++, p++, w >>= 1)
				{
					long lWavelet;
					
					// Given the currenet level l and subscript uPos, compute the sum of the portion in wavelet after uPos
					size_t uPosInWavelet = uPos % w;
					if( 0 == l )
						lWavelet = uPosInWavelet;
					else
					{
						if( uPosInWavelet < w / 2)
							lWavelet = (long)uPosInWavelet + 1;
						else
							lWavelet = (long)(w - 1 - uPosInWavelet);
					}
					vlWavelets.push_back( lWavelet );
				}
			}
			
			// now find the combination of the coefficients of all dimensions
			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				long lWavelet = 1;

				size_t uCoefId = 0;
				for(size_t d = 0, uBase = 0, uCoefBase = 1;
					d < UGetNrOfDims(); 
					uBase += vuDimMaxLevels[d], uCoefBase *= vuCoefLengths[d], d++, p++)
				{
					size_t uLevel = vuCoefDim2Level[p];
					/*
					// skip the level if it is larger than the threshold 
					if( uLevel >= vuDimMaxLevels[d] )
						continue;
					*/
					
					lWavelet *= vlWavelets[uBase + uLevel];
					uCoefId += vvuSubLevel2Coef[d][(this->vuDataDimLengths[d] - 1) * vuDimLevels[d] + uLevel];
				}

				// update the corresponding wavelet coeffcients
				for(size_t b = 0; b < this->UGetNrOfBins(); b++)
				{
					double dWeight = vdBinWeights[b];
					double dWavelet = dWeight * (double)lWavelet;
					double dOriginalCoef = this->DGetBinCoef(b, uCoefId);
					this->_SetBinCoef(b, uCoefId, dOriginalCoef + dWavelet);
				}
			}

			// update the DC component
			size_t uVolSize = 1;
			for(size_t d = 0; d < this->UGetNrOfDims(); d++)
				uVolSize *= this->vuDimLengths[d];

			for(size_t b = 0; b < this->UGetNrOfBins(); b++)
			{
				double dWeight = vdBinWeights[b];
				double dWavelet = dWeight * (double)uVolSize;
				double dOriginalCoef = this->DGetBinCoef(b, 0);
				this->_SetBinCoef(b, 0, dOriginalCoef - dWavelet);
			}
			#endif	// #if	WITH_BOUNDARY_AWARE_DWT	
			// ADD-BY-LEETEN 10/10/2012-BEGIN

			if( !bIsFinalizedWithoutWavelet )	
			{
				for(size_t b = 0; b < UGetNrOfBins(); b++)
				{	
				// ADD-BY-LEETEN 10/08/2012-END
					// ADD-BY-LEETEN 10/29/2012-BEGIN
					this->vcBinCoefs[b]._Finalize(dWaveletDenomiator);

					// ADD-BY-LEETEN 10/29/2012-END
				}	// ADD-BY-LEETEN 10/08/2012
			}	
			// ADD-BY-LEETEN 11/11/2012-BEGIN
			#else	// #if	!WITH_COEF_POOL
			{
				for(size_t c = 0; c < uNrOfUpdatingCoefs; c++)
														{
					vector<size_t> vuPoolSubs;
					_ConvertIndexToSub(c, vuPoolSubs, this->vuDimMaxLevels);

					// decide the wavelet weight
					double dWavelet = +1.0;
					// ADD-BY-LEETEN 12/16/2012-BEGIN
					if( !bIsFinalizedWithoutWavelet )	
					{
					// ADD-BY-LEETEN 12/16/2012-END
					for(size_t d = 0; d < vuPoolSubs.size(); d++)
					{
						size_t uLevel = vuPoolSubs[d];
						if( uLevel >= 1 )
							dWavelet *= (double)(1 << (uLevel - 1));
							// Wavelet *= (T)sqrt((double)(1 << (uLevel - 1) ));
					}

					dWavelet = sqrt(dWavelet);
					}
					double dWeight = ( !bIsFinalizedWithoutWavelet )?(dWavelet / dWaveletDenomiator):1.0;
					this->vcCoefPools[c]._Finalize( dWeight );
					// this->vcCoefPools[c]._Weight( Wavelet / dWaveletDenomiator );
				}
			}	
			#endif	// #if	!WITH_COEF_POOL
			// ADD-BY-LEETEN 11/11/2012-END

			_ShowMemoryUsage(false);

			// ADD-BY-LEETEN 12/15/2012-BEGIN
			uNrOfNonZeroCoefs = 0;
			uCountInFullArray = 0;
			uCountInSparseArray = 0;
			#if	!WITH_COEF_POOL	// ADD-BY-LEETEN 11/11/2012
			for(size_t b = 0; b < UGetNrOfBins(); b++)
			  {
			    size_t uF, uS;
			    this->vcBinCoefs[b]._GetArraySize(uF, uS, dWaveletThreshold);
			    uCountInFullArray += uF;
			    uCountInSparseArray += uS;
			  }
			// ADD-BY-LEETEN 11/11/2012-BEGIN
			#else	// #if	!WITH_COEF_POOL
			for(size_t c = 0; c < this->uNrOfUpdatingCoefs; c++)
			{
				size_t uF, uS;
				this->vcCoefPools[c]._GetArraySize(uF, uS, dWaveletThreshold);
				uCountInFullArray += uF;
				uCountInSparseArray += uS;
			}
			#endif	// #if	!WITH_COEF_POOL	
			uNrOfNonZeroCoefs = uCountInFullArray + uCountInSparseArray;
			// ADD-BY-LEETEN 12/15/2012-END
		}

		//! Compute and display statistics for the computed wavelet coefficients.
		virtual	
		void
		_ShowStatistics
		(
			void *_Reserved = NULL
		)
		{
			LOG_VAR(uCountInFullArray);
			LOG_VAR(uCountInSparseArray);
			LOG_VAR(uNrOfNonZeroCoefs);

			size_t uNrOfDataItems = this->uDataSize;

			double dCR = (double)(uNrOfDataItems * UGetNrOfBins()) / (double)uNrOfNonZeroCoefs;
			LOG_VAR(dCR);

			double dOverhead = (double)uNrOfNonZeroCoefs / (double)uNrOfDataItems;
			LOG_VAR(dOverhead);
		}

		//! Allocate the space to store coefficients for all bins. 
		/*! 
		*/
		virtual	// ADD-BY-LEETEN 09/29/2012
		void 
		_Allocate
		(
			void *_Reserved = NULL
		)
		{
			#if	!WITH_COEF_POOL	// ADD-BY-LEETEN 11/11/2012
			// ADD-BY-LEETEN 10/29/2012-BEGIN
			uNrOfCoefsInFullArray = min(
				(size_t)floor( (double)uSizeOfFullArrays/(double)(UGetNrOfBins() * sizeof(T))), 
				uNrOfCoefs);
			// ADD-BY-LEETEN 10/29/2012-END
			// ADD-BY-LEETEN 10/31/2012-BEGIN
			#if	WITH_1D_DIVISION
			// multiplied by the #coefficients s.t. later the indices can be computed without the extra multiplications
			for(size_t uBase = 1, d = 0; d < UGetNrOfDims(); uBase *= vuCoefLengths[d], d++)
				for(size_t c = 0; c < this->vvuSubLevel2Coef[d].size(); c++)
					this->vvuSubLevel2Coef[d][c] *= uBase;
			#endif	// #if	WITH_1D_DIVISION
			// ADD-BY-LEETEN 10/31/2012-END

			LOG_VAR(uSizeOfFullArrays);
			LOG_VAR(uNrOfCoefsInFullArray);

			// ADD-BY-LEETEN 10/29/2012-BEGIN
			vcBinCoefs.resize(uNrOfBins);
			for(size_t b = 0; b < uNrOfBins; b++)	
				vcBinCoefs[b]._Set(vuCoefLengths, uNrOfCoefsInFullArray);
			// ADD-BY-LEETEN 10/29/2012-END
			// ADD-BY-LEETEN 11/11/2012-BEGIN
			#else	// #if	!WITH_COEF_POOL
			// compute the #coef in the full array per bin SAT
			uNrOfCoefsInFullArray = min(
				(size_t)floor( (double)uSizeOfFullArrays/(double)(UGetNrOfBins() * sizeof(T))), 
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

			///////////////////////////////
			vector<size_t> vuPoolDimLengths;
			vuPoolDimLengths.resize(this->UGetNrOfDims());
			this->vcCoefPools.resize(this->uNrOfUpdatingCoefs);
			for(size_t c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				vector<size_t> vuPoolSubs;
				_ConvertIndexToSub(c, vuPoolSubs, this->vuDimLevels);

				bool bIsSparse = false;

				// decide the pool size
				for(size_t d = 0; d < this->UGetNrOfDims(); d++)
				{
					vuPoolDimLengths[d] = (!vuPoolSubs[d])?1:(1<<vuPoolSubs[d] - 1);
					// decide whether the array is sparse.
					if( vuPoolSubs[d] > vuOptimalDimLevel[d] )
						bIsSparse = true;
				}

				// ADD-BY-LEETEN 11/11/2012-BEGIN
				if( !uNrOfCoefsInFullArray )
				  bIsSparse = true;
				// ADD-BY-LEETEN 11/11/2012-END

				this->vcCoefPools[c]._Set(
					this->UGetNrOfBins(),
					vuPoolDimLengths,
					this->vuMaxCounts[c],	// ADD-BY-LEETEN 11/11/2012
					bIsSparse);
			}
			#endif	// #if	!WITH_COEF_POOL	
			// ADD-BY-LEETEN 11/11/2012-END

			_ShowMemoryUsage(false); // ADD-BY-LEETEN 11/14/2012
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
			#if	!WITH_COEF_POOL	// ADD-BY-LEETEN 11/11/2012
			// ADD-BY-LEETEN 10/29/2012-BEGIN
			#if !WITH_1D_DIVISION			
			vector<size_t> vuCoefPos;
			vuCoefPos.resize(UGetNrOfDims());
			#endif	// #if !WITH_1D_DIVISION	
			// ADD-BY-LEETEN 10/29/2012-END

			// for each bin, apply wavelet transform
			vdSums.clear();

			// ADD-BY-LEETEN 10/06/2012-BEGIN
			size_t uNrOfWavelets = 0;
			for(vector<size_t>::iterator 
				ivuDimMaxLevels = vuDimMaxLevels.begin();
				ivuDimMaxLevels != vuDimMaxLevels.end();
				ivuDimMaxLevels ++)
				uNrOfWavelets += *ivuDimMaxLevels;

			#if	!WITH_PRECOMPUTED_WAVELET_BASIS	// ADD-BY-LEETEN 10/18/2012
			vector<double> vdWaveletBasis;
			vdWaveletBasis.resize( uNrOfWavelets );
			// ADD-BY-LEETEN 10/18/2012-BEGIN
			#else	// #if	!WITH_PRECOMPUTED_WAVELET_BASIS
			vector<char> vbWaveletSigns;
			vbWaveletSigns.resize( uNrOfWavelets );
			#endif	// #if	!WITH_PRECOMPUTED_WAVELET_BASIS
			// ADD-BY-LEETEN 10/18/2012-END

			// ADD-BY-LEETEN 10/06/2012-END

			// ADD-BY-LEETEN 10/19/2012-BEGIN
			vector<size_t> vuPosLevelProduct;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				vuPosLevelProduct.push_back(vuPos[d] * vuDimLevels[d]);
			// ADD-BY-LEETEN 10/19/2012-END

			// ADD-BY-LEETEN 10/12/2012-BEGIN
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
					#if	!WITH_PRECOMPUTED_WAVELET_BASIS	// ADD-BY-LEETEN 10/18/2012
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
					// ADD-BY-LEETEN 10/18/2012-BEGIN
					#else	// #if	!WITH_PRECOMPUTED_WAVELET_BASIS	
					vbWaveletSigns[p] = ( uPos % w < w / 2 )?+1:-1;
					#endif	// #if	!WITH_PRECOMPUTED_WAVELET_BASIS	
					// ADD-BY-LEETEN 10/18/2012-END
				}
			}
			// ADD-BY-LEETEN 10/12/2012-END
		
			vdSums.resize(UGetNrOfBins());

			// now find the combination of the coefficients of all dimensions 
			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				#if WITH_1D_DIVISION
				size_t uCoefId = 0;
				#endif	// #if !WITH_1D_DIVISION
				double dWavelet = 1.0;
				#if	WITH_PRECOMPUTED_WAVELET_BASIS	
				int iWaveletSign = 1;
				#endif	// #if	WITH_PRECOMPUTED_WAVELET_BASIS	
				#if	!WITH_COEF_DIM_2_WAVELET	// ADD-BY-LEETEN 10/21/2012
				for(size_t d = 0, uBase = 0;
					d < UGetNrOfDims(); 
					uBase += vuDimMaxLevels[d], d++, p++)
				// ADD-BY-LEETEN 10/21/2012-BEGIN
				#else	// #if	!WITH_COEF_DIM_2_WAVELET	
				for(size_t d = 0; d < UGetNrOfDims(); d++, p++)
				#endif	// #if	!WITH_COEF_DIM_2_WAVELET	
				// ADD-BY-LEETEN 10/21/2012-END
				{
					#if WITH_1D_DIVISION	
					// update the index of the current coefficients
					uCoefId += vvuSubLevel2Coef[d][vuPosLevelProduct[d] + vuCoefDim2Level[p]];
					// ADD-BY-LEETEN 10/29/2012-BEGIN
					#else	// #if WITH_1D_DIVISION	
					vuCoefPos[d] = vvuSubLevel2Coef[d][vuPosLevelProduct[d] + vuCoefDim2Level[p]];
					#endif	// #if WITH_1D_DIVISION	
					// ADD-BY-LEETEN 10/29/2012-END

					// combine the wavelet basis value
					#if	!WITH_COEF_DIM_2_WAVELET	// ADD-BY-LEETEN 10/21/2012
					#if	!WITH_PRECOMPUTED_WAVELET_BASIS	
					dWavelet *= vdWaveletBasis[uBase + vuCoefDim2Level[p]];	
					#else	// #if	!WITH_PRECOMPUTED_WAVELET_BASIS	
					if( vbWaveletSigns[uBase + vuCoefDim2Level[p]] < 0 )
						iWaveletSign = -iWaveletSign;
					#endif	// #if	!WITH_PRECOMPUTED_WAVELET_BASIS	
					// ADD-BY-LEETEN 10/21/2012-BEGIN
					#else	// #if	!WITH_COEF_DIM_2_WAVELET	
					size_t uWavelet = vuCoefDim2Wavelet[p];
					#if	!WITH_PRECOMPUTED_WAVELET_BASIS	
					dWavelet *= vdWaveletBasis[uWavelet];
					#else	// #if	!WITH_PRECOMPUTED_WAVELET_BASIS	
					if( vbWaveletSigns[uWavelet] < 0 )
						iWaveletSign = -iWaveletSign;
					#endif	// #if	!WITH_PRECOMPUTED_WAVELET_BASIS	
					#endif	// #if	!WITH_COEF_DIM_2_WAVELET	
					// ADD-BY-LEETEN 10/21/2012-END
				}
				#if	WITH_PRECOMPUTED_WAVELET_BASIS	
				dWavelet = vdWaveletBasisPerUpdatingCoefs[c];
				if( iWaveletSign < 0 )
					dWavelet = -dWavelet;
				#endif	// #if	WITH_PRECOMPUTED_WAVELET_BASIS	

				for(size_t b = 0; b < UGetNrOfBins(); b++)
				{
					// ADD-BY-LEETEN 10/29/2012-BEGIN
					double dWaveletCoef;
					// ADD-BY-LEETEN 10/31/2012-BEGIN
					#if	WITH_1D_DIVISION	
					this->vcBinCoefs[b]._GetAtIndex(uCoefId, dWaveletCoef);
					#else	// #if	WITH_1D_DIVISION
					// ADD-BY-LEETEN 10/31/2012-END
					this->vcBinCoefs[b]._GetAtPos(vuCoefPos, dWaveletCoef);
					#endif	// #if	WITH_1D_DIVISION	// ADD-BY-LEETEN 10/31/2012
					// ADD-BY-LEETEN 10/29/2012-END

					if( fabs(dWaveletCoef) >= dWaveletThreshold )
					{
						// update the corresponding wavelet coeffcients
						double dIncremental = dWaveletCoef * dWavelet;
						vdSums[b] += dIncremental;
					}
				}
			}
			// ADD-BY-LEETEN 11/11/2012-BEGIN
			#else	// #if	!WITH_COEF_POOL	
			vector<size_t> vuEmpty;	// this empty vector is used to speed up the access of elements

			vdSums.resize(UGetNrOfBins());
			vector<double> vdWaveletBasis;
			vdWaveletBasis.resize( this->uNrOfWaveletsFromAllDims );

			vector<size_t> vuSubs;
			vuSubs.resize( this->uNrOfWaveletsFromAllDims );

			vector<size_t> vuPoolSubs;
			vuPoolSubs.resize( UGetNrOfDims() );

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
		
			// now find the combination of the coefficients of all dimensions 
			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				double dWavelet = 1.0;
				for(size_t d = 0, uBase = 0;
					d < UGetNrOfDims(); 
					uBase += vuDimMaxLevels[d], d++, p++)
				{
					vuPoolSubs[d] = vuSubs[uBase + vuCoefDim2Level[p]];
					dWavelet *= vdWaveletBasis[uBase + vuCoefDim2Level[p]];	
				}

				if( this->vcCoefPools[c].BIsSparse() )
				{
					vector< pair<size_t, double> > vpairCoefs;
					this->vcCoefPools[c]._GetCoefSparse
					(
						vuPoolSubs,
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
					size_t uIndex;
					for(size_t b = 0; b < UGetNrOfBins(); b++)
					{
						double dWaveletCoef;
						this->vcCoefPools[c]._GetAt(
							b, 
							( !b )?vuPoolSubs:vuEmpty, 
							uIndex, dWaveletCoef);

						if( fabs(dWaveletCoef) >= dWaveletThreshold )
							// update the corresponding wavelet coeffcients
							vdSums[b] += dWaveletCoef * dWavelet;
					}
				}
			}
			#endif	// #if	!WITH_COEF_POOL	
			// ADD-BY-LEETEN 11/11/2012-END
			for(size_t b = 0; b < UGetNrOfBins(); b++)
				vdSums[b] /= dWaveletDenomiator;
		}

		CWaveletSATEncoder():
			bIsFinalizedWithoutWavelet(false)
		{
		}
	};
}
