#pragma once

// MOD-BY-LEETEN 10/29/2012-FROM:	#define WITH_VECTORS_FOR_COUNTED_COEFS	1
#define WITH_VECTORS_FOR_COUNTED_COEFS	0
// MOD-BY-LEETEN 10/29/2012-END

// ADD-BY-LEETEN 10/10/2012-BEGIN
#define	WITH_BOUNDARY_AWARE_DWT		0
// ADD-BY-LEETEN 10/10/2012-END

// ADD-BY-LEETEN 10/29/2012-BEGIN
//! Decide whether the coefficients are stored in the class CSepDWTData
#define	WITH_SEP_DWT_DATA_CLASS		1
// ADD-BY-LEETEN 10/29/2012-END

#include <map>	
#if	0	// DEL-BY-LEETEN 10/29/2012-BEGIN
#if defined (WIN32)
	#include <psapi.h>	
	#pragma comment (lib, "psapi.lib")
#else	// #if defined (WIN32)
	#include <sys/time.h>
	#include <sys/resource.h>
#endif	// #if defined (WIN32)
#endif		// DEL-BY-LEETEN 10/29/2012-END

#include <vector>
using namespace std;
#include <math.h>

#include "utils.h"		// ADD-BY-LEETEN 10/29/2012
#include "SepDWTHeader.h"	
#include "SepDWTData.h"		// ADD-BY-LEETEN 10/29/2012
#include "EncoderBase.h"	

#include "liblog.h"	

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
// MOD-BY-LEETEN 10/25/2012-FROM:	namespace WaveletSAT
namespace WaveletSAT
// MOD-BY-LEETEN 10/25/2012-END
{
	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

	//! The class of our WaveletSAT algorithm.
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
		public CSepDWTHeader,
		public CEncoderBase<T, double>	
	{
protected:	
		#if	WITH_VECTORS_FOR_COUNTED_COEFS
		struct CTempCoef
		{
			size_t uCount;
			double dCoef;
			CTempCoef ():
				uCount(0),
				dCoef(0.0)
			{
			}
			
			CTempCoef(size_t uC, double dC):
				uCount(uC),
				dCoef(dC)
			{
			}
		};
		#endif	// #if	WITH_VECTORS_FOR_COUNTED_COEFS

		//! The flag whether the method _Finalize() should multiple the result by the wavelet coefficients
		bool bIsFinalizedWithoutWavelet;

		//! #Coefs stored in full arrays
		size_t uNrOfCoefsInFullArray;	

		#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012

		//! pool of the coefficents
		/*!
		*/
		vector< vector<double> > vvdBinCoefs;
		vector< map<size_t, double> > vmapBinCoefs;

		#if	WITH_VECTORS_FOR_COUNTED_COEFS
		vector< map<size_t, CTempCoef> > vmapBinTempCoefs;
		vector< vector< pair<size_t, double> > > vvdBinTempCoefs; 
		#endif	// #if	WITH_VECTORS_FOR_COUNTED_COEFS

		// ADD-BY-LEETEN 10/29/2012-BEGIN
		#else	// #if	!WITH_SEP_DWT_DATA_CLASS	
		vector< CSepDWTData<double> > vcBinCoefs;
		#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	
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
			// ADD-BY-LEETEN 10/29/2012-BEGIN
			#if	WITH_SEP_DWT_DATA_CLASS	
			vector<size_t> vuCoefPos;
			vuCoefPos.resize(UGetNrOfDims());
			#endif	// #if	WITH_SEP_DWT_DATA_CLASS	
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
					// MOD-BY-LEETEN 10/19/2012-FROM:	vlWavelets.push_back( lWavelet );
					vlWavelets[p] = lWavelet;
					// MOD-BY-LEETEN 10/19/2012-END
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
							
				#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012
				size_t uCoefId = 0;
				#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012

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

					#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012
					uCoefId += vvuSubLevel2Coef[d][vuPosLevelProduct[d] + uLevel];

					// ADD-BY-LEETEN 10/29/2012-BEGIN
					#else	// #if	!WITH_SEP_DWT_DATA_CLASS	
					vuCoefPos[d] = vvuSubLevel2Coef[d][vuPosLevelProduct[d] + uLevel];
					#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	
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
				#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012
				size_t uCoefId = 0;
				#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012
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
					#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012-BEGIN
					uCoefId += vvuSubLevel2Coef[d][vuPosLevelProduct[d] + uLevel];
					
					// ADD-BY-LEETEN 10/29/2012-BEGIN
					#else	// #if	!WITH_SEP_DWT_DATA_CLASS	
					vuCoefPos[d] = vvuSubLevel2Coef[d][vuPosLevelProduct[d] + uLevel]
					#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	
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

				#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012
				// update the corresponding wavelet coeffcients
				if( uCoefId < uNrOfCoefsInFullArray ) 
				{
					vvdBinCoefs[uBin][uCoefId] += dWavelet;
				}
				else 
				{
					#if	!WITH_VECTORS_FOR_COUNTED_COEFS	
					map<size_t, double>::iterator ipairCoef = this->vmapBinCoefs[uBin].find(uCoefId);
					if(this->vmapBinCoefs[uBin].end() == ipairCoef )
					{
						this->vmapBinCoefs[uBin].insert(pair<size_t, double>(uCoefId, dWavelet));
					#else	// #if	!WITH_VECTORS_FOR_COUNTED_COEFS
#if 0 // MOD-BY-LEETEN 10/26/2012-FROM:
					map<size_t, CTempCoef>::iterator ipairTempCoef = this->vmapBinTempCoefs[uBin].find(uCoefId);
#else // MOD-BY-LEETEN 10/26/2012-TO:
					typename map<size_t, CTempCoef>::iterator ipairTempCoef = this->vmapBinTempCoefs[uBin].find(uCoefId);
#endif // MOD-BY-LEETEN 10/26/2012-END
					if(this->vmapBinTempCoefs[uBin].end() == ipairTempCoef )
					{
						this->vmapBinTempCoefs[uBin].insert(pair<size_t, CTempCoef>(uCoefId, CTempCoef(1, dWavelet)));
					#endif	// #if	!WITH_VECTORS_FOR_COUNTED_COEFS

						static size_t uCount;
						static size_t uMaxCount = 100000;
						if( 0 == uCount % uMaxCount )
						{
							LOG_VAR(uCount);
							#if defined(WIN32)
							PROCESS_MEMORY_COUNTERS memCounter;
							BOOL result = GetProcessMemoryInfo(
									GetCurrentProcess(),
									&memCounter,
									sizeof( memCounter ));
							LOG_VAR(memCounter.WorkingSetSize);
							#else	// #if defined(WIN32)
							int who = RUSAGE_SELF; 
							struct rusage usage; 
							int ret; 
							getrusage(who,&usage);
							LOG_VAR(usage.ru_maxrss);
							#endif	// #if defined(WIN32)
						}
						uCount++;
					}
					else
					{
						#if	!WITH_VECTORS_FOR_COUNTED_COEFS	
						ipairCoef->second += dWavelet;
						#else	// #if	!WITH_VECTORS_FOR_COUNTED_COEFS
						ipairTempCoef->second.uCount++;
						ipairTempCoef->second.dCoef += dWavelet;

						if( this->vuMaxCounts[c] == ipairTempCoef->second.uCount )
						{
							this->vvdBinTempCoefs[uBin].push_back(pair<size_t, double>(uCoefId, ipairTempCoef->second.dCoef));
							this->vmapBinTempCoefs[uBin].erase(ipairTempCoef);
						}
						#endif	// #if	!WITH_VECTORS_FOR_COUNTED_COEFS
					}
				}
				// ADD-BY-LEETEN 10/29/2012-BEGIN
				#else	// #if	!WITH_SEP_DWT_DATA_CLASS	
				this->vcBinCoefs[uBin]._AddAtPos(vuCoefPos, dWavelet);
				#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	
				// ADD-BY-LEETEN 10/29/2012-END
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
			CHeaderBase::_SetLong(eName, lValue);
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

		#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012
		// ADD-BY-LEETEN 10/08/2012-BEGIN
		virtual 
		double 
		DGetBinCoef
		(
			size_t uBin,
			size_t uCoefId,
			void *_Reserved = NULL
		)
		{
			double dCoef = 0.0;
			if( uCoefId < this->uNrOfCoefsInFullArray )
				dCoef = vvdBinCoefs[uBin][uCoefId];
			else
			{
				map<size_t, double>::iterator ipairCoef = this->vmapBinCoefs[uBin].find(uCoefId);
				if( this->vmapBinCoefs[uBin].end() != ipairCoef )
					dCoef = ipairCoef->second;
			}
			return dCoef;
		}

		virtual 
		void
		_SetBinCoef
		(
			size_t uBin,
			size_t uCoefId,
			double dCoef, 
			void *_Reserved = NULL
		)
		{
			if( uCoefId < this->uNrOfCoefsInFullArray )
				vvdBinCoefs[uBin][uCoefId] = dCoef;
			else
			{
				map<size_t, double>::iterator ipairCoef = this->vmapBinCoefs[uBin].find(uCoefId);
				if( this->vmapBinCoefs[uBin].end() != ipairCoef )
				// ADD-BY-LEETEN 10/12/2012-BEGIN
				{	
					if( !dCoef )
						this->vmapBinCoefs[uBin].erase(ipairCoef);
					else
				// ADD-BY-LEETEN 10/12/2012-END
					ipairCoef->second = dCoef;
				}	// ADD-BY-LEETEN 10/12/2012
				else
				if( dCoef )	// ADD-BY-LEETEN 10/10/2012
					this->vmapBinCoefs[uBin].insert(pair<size_t, double>(uCoefId, dCoef));
			}
		}
		// ADD-BY-LEETEN 10/08/2012-END
		#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012

		//! Finalize the computation of SAT
		virtual	
		void 
		_Finalize
		(
			void *_Reserved = NULL
		)
		{
			_ShowMemoryUsage();	// ADD-BY-LEETEN 10/29/2012

			vector<size_t> vuSub;	// ADD-BY-LEETEN 10/06/2012

			#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012
			for(size_t b = 0; b < UGetNrOfBins(); b++)
			{	
				// ADD-BY-LEETEN 10/08/2012-BEGIN
				#if	WITH_VECTORS_FOR_COUNTED_COEFS
#if 0 // MOD-BY-LEETEN 10/26/2012-FROM:
				for(map<size_t, CTempCoef>::iterator 
#else // MOD-BY-LEETEN 10/26/2012-TO:
				for(typename map<size_t, CTempCoef>::iterator 
#endif // MOD-BY-LEETEN 10/26/2012-END
					ipairTempCoef = this->vmapBinTempCoefs[b].begin();
					ipairTempCoef != this->vmapBinTempCoefs[b].end();
					ipairTempCoef++)
					this->vmapBinCoefs[b].insert(pair<size_t, double>(ipairTempCoef->first, ipairTempCoef->second.dCoef));
				this->vmapBinTempCoefs[b].clear();

				for(vector< pair<size_t, double> >::iterator 
					ivdTempCoef = this->vvdBinTempCoefs[b].begin();
					ivdTempCoef != this->vvdBinTempCoefs[b].end();
					ivdTempCoef++)
					this->vmapBinCoefs[b].insert(*ivdTempCoef);
				this->vvdBinTempCoefs[b].clear();
				#endif	// #if	!WITH_VECTORS_FOR_COUNTED_COEFS
			}	// ADD-BY-LEETEN 10/10/2012
			#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012

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
					#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012
					for(size_t w = 0; w < this->vvdBinCoefs[b].size(); w++)
					{
						double dCoef = this->vvdBinCoefs[b][w];
						if( dCoef )
						{
							double dWavelet = +1.0;
							_ConvetIndexToSub(w, vuSub, vuCoefLengths);

							for(size_t d = 0; d < vuSub.size(); d++)
							{
								size_t uSub = vuSub[d];
								if( uSub >= 1 )
								{
									size_t uLevel = (size_t)ceil(log( (double)(uSub + 1) ) / log(2.0) );
									dWavelet *= sqrt((double)(1 << (uLevel - 1) ));
								}
							}
							this->vvdBinCoefs[b][w] *= dWavelet / dWaveletDenomiator;
						}
					}

					for(map<size_t, double>::iterator 
						ipairCoef = this->vmapBinCoefs[b].begin();
						ipairCoef != this->vmapBinCoefs[b].end();
						ipairCoef++)
					{
						double dCoef = ipairCoef->second;
						if( dCoef )
						{
							double dWavelet = 1.0;

							_ConvetIndexToSub(ipairCoef->first, vuSub, vuCoefLengths);

							for(size_t d = 0; d < vuSub.size(); d++)
							{
								size_t uSub = vuSub[d];
								if( uSub >= 1 )
								{
									size_t uLevel = (size_t)ceil(log( (double)(uSub + 1) ) / log(2.0) );
									dWavelet *= sqrt((double)(1 << (uLevel - 1) ));	
								}
							}
							ipairCoef->second *= dWavelet / dWaveletDenomiator;
						}
					}
					// ADD-BY-LEETEN 10/29/2012-BEGIN
					#else	// #if	!WITH_SEP_DWT_DATA_CLASS	
					this->vcBinCoefs[b]._Finalize(dWaveletDenomiator);

					#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	
					// ADD-BY-LEETEN 10/29/2012-END
				}	// ADD-BY-LEETEN 10/08/2012
			}	
		}

		//! Compute and display statistics for the computed wavelet coefficients.
		virtual	
		void
		_ShowStatistics
		(
			void *_Reserved = NULL
		)
		{
			size_t uNrOfNonZeroCoefs = 0;
			#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012
			for(size_t b = 0; b < UGetNrOfBins(); b++)
			{
				double dEnergy = 0.0;
				for(size_t w = 0; w < this->vvdBinCoefs[b].size(); w++)
				{
					double dCoef = this->vvdBinCoefs[b][w];
					dEnergy += pow(dCoef, 2.0);
					if( fabs(dCoef) > this->dWaveletThreshold )
						uNrOfNonZeroCoefs++;
				}
				for(map<size_t, double>::iterator
					ipairCoef = vmapBinCoefs[b].begin();
					ipairCoef != vmapBinCoefs[b].end();
					ipairCoef++)
				{
					double dCoef = ipairCoef->second;
					dEnergy += pow(dCoef, 2.0);
					if( fabs(dCoef) > this->dWaveletThreshold )
						uNrOfNonZeroCoefs++;
				}
				// printf("Energy[%d] = %f\n", b, dEnergy);
			}
			// ADD-BY-LEETEN 10/29/2012-BEGIN
			#else	// #if	!WITH_SEP_DWT_DATA_CLASS	
			for(size_t b = 0; b < UGetNrOfBins(); b++)
			{
				size_t uCount;
				this->vcBinCoefs[b]._GetNrOfNonZeroCoefs(uCount, dWaveletThreshold);
				uNrOfNonZeroCoefs += uCount;
			}
			#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	
			// ADD-BY-LEETEN 10/29/2012-END

			LOG_VAR(uNrOfNonZeroCoefs);

			size_t uNrOfDataItems = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				LOG_VAR(vuDimLengths[d]);
				uNrOfDataItems *= CHeaderBase::vuDimLengths[d];
			}
			LOG_VAR(uNrOfDataItems);
			LOG_VAR(UGetNrOfBins());

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
			// ADD-BY-LEETEN 10/29/2012-BEGIN
			#if	!WITH_SEP_DWT_DATA_CLASS	
			// multiplied by the #coefficients s.t. later the indices can be computed without the extra multiplications
			for(size_t uBase = 1, d = 0; d < UGetNrOfDims(); uBase *= vuCoefLengths[d], d++)
				for(size_t c = 0; c < this->vvuSubLevel2Coef[d].size(); c++)
					this->vvuSubLevel2Coef[d][c] *= uBase;
			// ADD-BY-LEETEN 10/29/2012-END

			uNrOfCoefsInFullArray = min(
				(size_t)floor( (double)uSizeOfFullArrays/(double)(UGetNrOfBins() * sizeof(this->vvdBinCoefs[0][0]))), 
				uNrOfCoefs);

			// ADD-BY-LEETEN 10/29/2012-BEGIN
			#else	// #if	!WITH_SEP_DWT_DATA_CLASS	
			uNrOfCoefsInFullArray = min(
				(size_t)floor( (double)uSizeOfFullArrays/(double)(UGetNrOfBins() * sizeof(T))), 
				uNrOfCoefs);
			#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	
			// ADD-BY-LEETEN 10/29/2012-END

			LOG_VAR(uSizeOfFullArrays);
			LOG_VAR(uNrOfCoefsInFullArray);

			#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012
			vvdBinCoefs.resize(uNrOfBins);
			for(size_t b = 0; b < uNrOfBins; b++)
				vvdBinCoefs[b].resize(uNrOfCoefsInFullArray);

			for(size_t b = 0; b < uNrOfBins; b++)
				vmapBinCoefs.resize(uNrOfBins);

			#if	WITH_VECTORS_FOR_COUNTED_COEFS
			for(size_t b = 0; b < uNrOfBins; b++)	
			{
				vmapBinTempCoefs.resize(uNrOfBins);
				vvdBinTempCoefs.resize(uNrOfBins);
			}
			#endif	// #if	WITH_VECTORS_FOR_COUNTED_COEFS

			// ADD-BY-LEETEN 10/29/2012-BEGIN
			#else	// #if	!WITH_SEP_DWT_DATA_CLASS	
			vcBinCoefs.resize(uNrOfBins);
			for(size_t b = 0; b < uNrOfBins; b++)	
				vcBinCoefs[b]._Set(vuCoefLengths, uNrOfCoefsInFullArray);
			
			#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	
			// ADD-BY-LEETEN 10/29/2012-END
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
			// ADD-BY-LEETEN 10/29/2012-BEGIN
			#if	WITH_SEP_DWT_DATA_CLASS	
			vector<size_t> vuCoefPos;
			vuCoefPos.resize(UGetNrOfDims());
			#endif	// #if	WITH_SEP_DWT_DATA_CLASS	
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
				#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012
				size_t uCoefId = 0;
				#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012
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
					#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012-BEGIN
					// update the index of the current coefficients
					uCoefId += vvuSubLevel2Coef[d][vuPosLevelProduct[d] + vuCoefDim2Level[p]];
					// ADD-BY-LEETEN 10/29/2012-BEGIN
					#else	// #if	!WITH_SEP_DWT_DATA_CLASS	
					vuCoefPos[d] = vvuSubLevel2Coef[d][vuPosLevelProduct[d] + vuCoefDim2Level[p]];
					#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	
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
					#if	!WITH_SEP_DWT_DATA_CLASS	// ADD-BY-LEETEN 10/29/2012
					double dWaveletCoef = DGetBinCoef(b, uCoefId);
					// ADD-BY-LEETEN 10/29/2012-BEGIN
					#else	// #if	!WITH_SEP_DWT_DATA_CLASS	
					double dWaveletCoef;
					this->vcBinCoefs[b]._GetAtPos(vuCoefPos, dWaveletCoef);
					#endif	// #if	!WITH_SEP_DWT_DATA_CLASS	
					// ADD-BY-LEETEN 10/29/2012-END

					if( fabs(dWaveletCoef) >= dWaveletThreshold )	// MOD-BY-LEETEN 10/06/2012-FROM:	if( 0.0 != dWaveletCoef )
					{
						// update the corresponding wavelet coeffcients
						double dIncremental = dWaveletCoef * dWavelet;
						vdSums[b] += dIncremental;
					}
				}
			}
			for(size_t b = 0; b < UGetNrOfBins(); b++)
				vdSums[b] /= dWaveletDenomiator;
		}

		// MOD-BY-LEETEN 10/25/2012-FROM:	CBase():
		CWaveletSATEncoder():
		// MOD-BY-LEETEN 10/25/2012-END
			bIsFinalizedWithoutWavelet(false)
		{
		}
	};
}
