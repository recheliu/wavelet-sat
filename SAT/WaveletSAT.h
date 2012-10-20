#pragma once

#define WITH_VECTORS_FOR_COUNTED_COEFS	1

#if	0	// MOD-BY-LEETEN 10/19/2012-FROM:
#define	WITH_LOG_WAVELET_BASIS		1	// ADD-BY-LEETEN 10/18/2012
#else		// MOD-BY-LEETEN 10/19/2012-TO:
#define	WITH_LOG_WAVELET_BASIS		0
#endif		// MOD-BY-LEETEN 10/19/2012-END

// ADD-BY-LEETEN 10/10/2012-BEGIN
#define	WITH_BOUNDARY_AWARE_DWT		0
// ADD-BY-LEETEN 10/10/2012-END

#include <map>	
#if defined (WIN32)
	#include <psapi.h>	
	#pragma comment (lib, "psapi.lib")
#else	// #if defined (WIN32)
	#include <sys/time.h>
	#include <sys/resource.h>
#endif	// #if defined (WIN32)

#include <vector>
using namespace std;
#include <math.h>

#include "Base.h"	
#include "liblog.h"	

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
  struct CWaveletTempCoef 
  {
    size_t uCount;
    double dCoef;
    CWaveletTempCoef ():
    uCount(0),
      dCoef(0.0)
    {}
    CWaveletTempCoef (size_t uC, double dC):
    uCount(uC),
      dCoef(dC)
    {}
  };

	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

	//! The base class of WaveletSAT.
	/*!
	In order to apply wavelet transform for SATs or Integral Histograms, please follow the procedure below:
	1. Setup the data size by _SetDimLengths().
	2. Setup the #SATs by _AllocateBins().
	3. Call _Update() for each data point to incrementally compute wavelet.
	4. Call _Finalize() to finish the computation.

	To query the entry for a certain location:
	1. Call _GetAllSums() to get all SAT values for the given location.
	*/
	template<class T>
	class CBase
		:public SAT::CBase<T>	// ADD-BY-LEETEN 09/29/2012
	{
protected:	
		// ADD-BY-LEETEN 10/19/2012-BEGIN
		//! Total #coefficients
		size_t uNrOfCoefs;
		// ADD-BY-LEETEN 10/19/2012-END

		//! The flag whether the method _Finalize() should multiple the result by the wavelet coefficients
		bool bIsFinalizedWithoutWavelet;

		//! The maximun number of wavelet levels per dim, which are specified by the user.
		/*!
		m[0], ... m[d], ... m[D - 1]
		*/
		vector<size_t> vuDimMaxLevels;
		
		//! #Coefs per dim.
		/*!
		c[0], ..., c[d], ... c[D - 1]
		*/
		vector<size_t> vuCoefLengths;

		//! pool of the coefficents
		/*!
		*/
		vector< vector<double> > vvdBinCoefs;
		vector< map<size_t, double> > vmapBinCoefs;
		//! #Coefs stored in full arrays
		size_t uNrOfCoefsInFullArray;	

		// ADD-BY-LEETEN 10/10/2012-BEGIN
		//! The weight sum for each bin
		vector<double> vdBinWeights;
		// ADD-BY-LEETEN 10/10/2012-END

		//! Total #coefs stored in full arrays from all bin SATs.
		size_t uSizeOfFullArrays;

		#if	WITH_VECTORS_FOR_COUNTED_COEFS
		vector< map<size_t, CWaveletTempCoef> > vmapBinTempCoefs;
		vector< vector< pair<size_t, double> > > vvdBinTempCoefs; 
		vector<size_t> vuMaxCounts;	// ADD-BY-LEETEN 10/19/2012
		#endif	// #if	WITH_VECTORS_FOR_COUNTED_COEFS
		
		//! The dim length of the input data
		vector<size_t> vuDataDimLengths;

		//! The dim length of the coefficients
		/*!
		n[0], ..., n[d], ... n[D - 1]
		*/
		vector<size_t> vuDimLengths;
		
		//! #levels per dim
		/*!
		l[0] = log2(n[0]) + 1, ... l[d] = log2(n[d]) + 1, ... l[D - 1] = log2(n[D - 1]) + 1
		*/
		vector<size_t> vuDimLevels;

		//! Product of influencing numbers of coefficients from all dimensions
		/*! 
		#coefficients to update ...
		C = prod l[d]
		*/
		size_t uNrOfUpdatingCoefs;

		//! #Dimensions
		/*!
		D
		*/
		// MOD-BY-LEETEN 10/17/2012-FROM:		size_t UGetNrOfDims() const {	return vuDimLengths.size();	};
		size_t 
		UGetNrOfDims
		(
			void *_Reserved = NULL
		) const 
		{	
			return vuDimLengths.size();	
		};
		// MOD-BY-LEETEN 10/17/2012-END
		
		size_t
		UConvetSubToIndex
		(
			const vector<size_t>& vuSub,
			void* _Reserved = NULL
		)
		{
			size_t uIndex = 0;
			for(size_t d = 0, vuSubSize = 1; d < vuSub.size(); vuSubSize *= this->vuDimLengths[d], d++)
				uIndex += vuSubSize * vuSub[d];
			return uIndex;
		}

		void
		_ConvetIndexToSub
		(
			size_t uIndex,
			vector<size_t>& vuSub,
			void* _Reserved = NULL
		)
		{
			if( this->UGetNrOfDims() != vuSub.size() )
				vuSub.resize(this->UGetNrOfDims());
			for(size_t d = 0; d < this->UGetNrOfDims(); d++)
			{
				vuSub[d] = uIndex % this->vuDimLengths[d];
				uIndex /= this->vuDimLengths[d];
			}
		}

		//! #Dimensions
		/*!
		B
		*/
		// MOD-BY-LEETEN 10/17/2012-FROM:	size_t UGetNrOfBins() const {	return vvdBinCoefs.size();	};
		size_t 
		UGetNrOfBins
		(
			void *_Reserved = NULL
		) const 
		{	
			return vvdBinCoefs.size();	
		};
		// MOD-BY-LEETEN 10/17/2012-END

		//! A lookup table to map the coefficient to its levels per dim.
		/*! size: D x C: 
		Map each coefficients c, c = 0 ... C - 1 = prod l[d] - 1, to the the corresponding D levels
		*/
		vector<size_t> vuCoefDim2Level;

		//! The look up table per dimension to quickly conver dim. subscript to the coefficient indices
		/*!
		For each dim d, 
			the size is l[d] x n[d]: map the spatial location i, i = 0 ... n[d], to its l[d] coefficient locations
		*/
		vector< vector<size_t> > vvuSubLevel2Coef;

		//! The denomator for Wavelet basis
		/*!
		W = sqrt(prod n[0], ... n[d], ... n[D])
		*/
		double dWaveletDenomiator;

		// ADD-BY-LEETEN 10/06/2012-BEGIN
		//! The mininal value of the projected wavelet coefficients
		double dWaveletThreshold;	
		// ADD-BY-LEETEN 10/06/2012-END

		//! Update the specified bin.
		void 
		_UpdateBin
		(
			const vector<size_t>& vuPos, 
			const T& value,
			size_t uBin, 
			double dWeight,
			void *_Reserved = NULL
		)
		{
			vdBinWeights[uBin] += dWeight;	// ADD-BY-LEETEN 10/10/2012

			#if	!WITH_LOG_WAVELET_BASIS	// ADD-BY-LEETEN 10/18/2012
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
			// ADD-BY-LEETEN 10/19/2012-BEGIN
			vector<size_t> vuPosLevelProduct;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				vuPosLevelProduct.push_back(vuPos[d] * vuDimLevels[d]);
			// ADD-BY-LEETEN 10/19/2012-END
			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				long lWavelet = 1;

				#if	WITH_VECTORS_FOR_COUNTED_COEFS
				// MOD-BY-LEETEN 10/19/2012-FROM:	size_t uMaxCountPerCoef = 1;	// max # coefficnet
				size_t uMaxCountPerCoef = this->vuMaxCounts[c];
				// MOD-BY-LEETEN 10/19/2012-END
				#endif	// #if	WITH_VECTORS_FOR_COUNTED_COEFS
				size_t uCoefId = 0;
				#if	0	// MOD-BY-LEETEN 10/19/2012-FROM:
				for(size_t d = 0, uBase = 0, uCoefBase = 1;
					d < UGetNrOfDims(); 
					uBase += vuDimMaxLevels[d], uCoefBase *= vuCoefLengths[d], d++, p++)
				#else		// MOD-BY-LEETEN 10/19/2012-TO:
				for(size_t d = 0, uBase = 0;
					d < UGetNrOfDims(); 
					uBase += vuDimMaxLevels[d], d++, p++)
				#endif		// MOD-BY-LEETEN 10/19/2012-END
				{
					size_t uLevel = vuCoefDim2Level[p];
					
					#if	0	// DEL-BY-LEETEN 10/19/2012-BEGIN
					#if	WITH_VECTORS_FOR_COUNTED_COEFS
					uMaxCountPerCoef *= (size_t)1 << (vuDimLevels[d] - uLevel);
					#endif	// #if	WITH_VECTORS_FOR_COUNTED_COEFS
					#endif		// DEL-BY-LEETEN 10/19/2012-END
					/*
					// skip the level if it is larger than the threshold 
					if( uLevel >= vuDimMaxLevels[d] )
						continue;
					*/
					
					// ADD-BY-LEETEN 10/18/2012-BEGIN
					if( 0 == d )
						lWavelet = vlWavelets[uBase + uLevel];
					else
					// ADD-BY-LEETEN 10/18/2012-END
					lWavelet *= vlWavelets[uBase + uLevel];

					#if	0	// MOD-BY-LEETEN 10/19/2012-FROM:
 					size_t uCoef = vvuSubLevel2Coef[d][vuPos[d] * vuDimLevels[d] + uLevel];
					uCoefId += uCoef * uCoefBase;
					#else		// MOD-BY-LEETEN 10/19/2012-TO:
					uCoefId += vvuSubLevel2Coef[d][vuPosLevelProduct[d] + uLevel];
					#endif		// MOD-BY-LEETEN 10/19/2012-END

					// ADD-BY-LEETEN 10/18/2012-BEGIN
					if( !lWavelet )
					{
						p += UGetNrOfDims() - d;
						break;
					}
					// ADD-BY-LEETEN 10/18/2012-END
				}
				// ADD-BY-LEETEN 10/18/2012-BEGIN
				if( !lWavelet )
					continue;
				// ADD-BY-LEETEN 10/18/2012-END
				double dWavelet = dWeight * (double)lWavelet;

			// ADD-BY-LEETEN 10/18/2012-BEGIN
			#else	// #if	!WITH_LOG_WAVELET_BASIS	
			vector<double> vdWaveletLogs;
			vector<char> vbWaveletSigns;
			// ADD-BY-LEETEN 10/19/2012-BEGIN
			vdWaveletLogs.resize(UGetNrOfDims() * vuDimLevels[0]);
			vbWaveletSigns.resize(UGetNrOfDims() * vuDimLevels[0]);
			// ADD-BY-LEETEN 10/19/2012-END

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
					char bSign;
					if( 0 == l )
					{
						lWavelet = (long)w / 2 - (long)uPosInWavelet;
						bSign = (lWavelet)?+1:0;
					}
					else
					{
						if( uPosInWavelet < w / 2)
							lWavelet = (long)uPosInWavelet;
						else
							lWavelet = (long)(w - uPosInWavelet);
						bSign = (lWavelet)?-1:0;
					}
					double dLogWavelet = (lWavelet)?log((double)lWavelet):0.0;
					#if	0	// MOD-BY-LEETEN 10/19/2012-FROM:
					vdWaveletLogs.push_back( dLogWavelet );
					vbWaveletSigns.push_back( bSign );
					#else		// MOD-BY-LEETEN 10/19/2012-TO:
					vdWaveletLogs[p] = dLogWavelet;
					vbWaveletSigns[p] = bSign;
					#endif		// MOD-BY-LEETEN 10/19/2012-END
				}
			}
			
			// now find the combination of the coefficients of all dimensions
			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				double dWaveletLog = 0.0;
				char bWaveletSign = 1;

				#if	WITH_VECTORS_FOR_COUNTED_COEFS
				size_t uMaxCountPerCoef = 1;	// max # coefficnet
				#endif	// #if	WITH_VECTORS_FOR_COUNTED_COEFS
				size_t uCoefId = 0;
				for(size_t d = 0, uBase = 0, uCoefBase = 1;
					d < UGetNrOfDims(); 
					uBase += vuDimMaxLevels[d], uCoefBase *= vuCoefLengths[d], d++, p++)
				{
					size_t uLevel = vuCoefDim2Level[p];
					#if	WITH_VECTORS_FOR_COUNTED_COEFS
					uMaxCountPerCoef *= (size_t)1 << (vuDimLevels[d] - uLevel);
					#endif	// #if	WITH_VECTORS_FOR_COUNTED_COEFS
					
					dWaveletLog += vdWaveletLogs[uBase + uLevel];
					bWaveletSign *= vbWaveletSigns[uBase + uLevel];
					if(!bWaveletSign)
					{
						p += UGetNrOfDims() - d;
						break;
					}

					size_t uCoef = vvuSubLevel2Coef[d][vuPos[d] * vuDimLevels[d] + uLevel];
					uCoefId += uCoef * uCoefBase;
				}
				if( !bWaveletSign )
					continue;

				double dWavelet = 0.0;
				dWavelet = dWeight * exp(dWaveletLog);
				if( bWaveletSign < 0)
					dWavelet = -dWavelet;
			#endif	// #if	!WITH_LOG_WAVELET_BASIS	
			// ADD-BY-LEETEN 10/18/2012-END

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
				map<size_t, CWaveletTempCoef>::iterator ipairTempCoef = this->vmapBinTempCoefs[uBin].find(uCoefId);
				if(this->vmapBinTempCoefs[uBin].end() == ipairTempCoef )
				{
					this->vmapBinTempCoefs[uBin].insert(pair<size_t, CWaveletTempCoef>(uCoefId, CWaveletTempCoef(1, dWavelet)));
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

					if( uMaxCountPerCoef == ipairTempCoef->second.uCount )
					{
						this->vvdBinTempCoefs[uBin].push_back(pair<size_t, double>(uCoefId, ipairTempCoef->second.dCoef));
						this->vmapBinTempCoefs[uBin].erase(ipairTempCoef);
					}
					#endif	// #if	!WITH_VECTORS_FOR_COUNTED_COEFS
				}
				}
			}
		}

		//! Given the subscript of a point, update the corresponding bin SAT
		/*!
		Note: As this method is not public, it is not intended to be directly used by the application. 
		Instead, the sub class should define additional methods to call this method.
		This design is for the case that not all data's area available
		*/
		virtual	// ADD-BY-LEETEN 09/29/2012
		void 
		_Update
		(
			const vector<size_t>& vuPos,
			const T& value,
			void *_Reserved = NULL
		)
		{
			vector< pair<size_t, double> > vpBins;

			_MapValueToBins(vuPos, value, vpBins);

			for(vector< pair<size_t, double> >::iterator ivpBin  = vpBins.begin(); ivpBin != vpBins.end(); ivpBin++)
				_UpdateBin(vuPos, value, ivpBin->first, ivpBin->second);
			
		}

		////////////////////////////////////////////////////////////////////
		/*
		The interface that should be overloaed by the sub class
		*/

		//! This method should be overloaded to return the indices of bins that will be changed
		virtual 
		void 
		_MapValueToBins
		(
			const vector<size_t>& vuPos,
			const T& value, 
			vector< pair<size_t, double> >& vuBins,
			void *_Reserved = NULL
		) = 0;
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			//! Total size in MB of all bin SATs.
			SIZE_OF_FULL_ARRAYS,

			//! Indicate whether the function _Finalize() should weight the result by wavelet coefficients
			FINALIZED_WITHOUT_WAVELET,

			NR_OF_ENUMS
		};

		virtual	
		void
		_SetLong(
			int eName,
			long lValue,
			void* _Reserved = NULL
		)
		{
			switch(eName)
			{
			case SIZE_OF_FULL_ARRAYS:
				uSizeOfFullArrays = (size_t)lValue * 1024 * 1024;
				break;
			}
		}

		virtual	
		void
		_SetBoolean(
			enum EParameter eName,
			bool bValue,
			void* _Reserved = NULL
		)
		{
			switch(eName)
			{
			case FINALIZED_WITHOUT_WAVELET:
				bIsFinalizedWithoutWavelet = bValue;
				break;
			}
		}

		//! Retunr the smallest value. 
		/*!
		This can be used to filter too small value caused by numerical error.
		*/
		double 
		DGetThreshold
		(
			void *_Reserved = NULL
		)
		{
			return dWaveletThreshold;
		}

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

		//! Finalize the computation of SAT
		virtual	
		void 
		_Finalize
		(
			void *_Reserved = NULL
		)
		{
			int iWaveletDenomiator = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				iWaveletDenomiator *= 1 << (vuDimLevels[d] - 1);
			dWaveletDenomiator = sqrt((double)iWaveletDenomiator);

			vector<size_t> vuSub;	// ADD-BY-LEETEN 10/06/2012

			// DEL-BY-LEETEN 10/08/2012:	if( !bIsFinalizedWithoutWavelet )	
			for(size_t b = 0; b < UGetNrOfBins(); b++)
			{	
				// ADD-BY-LEETEN 10/08/2012-BEGIN
				#if	WITH_VECTORS_FOR_COUNTED_COEFS
				for(map<size_t, CWaveletTempCoef>::iterator 
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
					#if	0	// MOD-BY-LEETEN 10/19/2012-FROM:
					size_t uCoef = vvuSubLevel2Coef[d][(this->vuDataDimLengths[d] - 1) * vuDimLevels[d] + uLevel];
					uCoefId += uCoef * uCoefBase;
					#else		// MOD-BY-LEETEN 10/19/2012-TO:
					uCoefId += vvuSubLevel2Coef[d][(this->vuDataDimLengths[d] - 1) * vuDimLevels[d] + uLevel];
					#endif		// MOD-BY-LEETEN 10/19/2012-END
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

			#if	0	// MOD-BY-LEETEN 10/10/2012-FROM:
				if( !bIsFinalizedWithoutWavelet )	
				{
			#else		// MOD-BY-LEETEN 10/10/2012-TO:	
			if( !bIsFinalizedWithoutWavelet )	
			{
				for(size_t b = 0; b < UGetNrOfBins(); b++)
				{	
			#endif		// MOD-BY-LEETEN 10/10/2012-END
				// ADD-BY-LEETEN 10/08/2012-END
				for(size_t w = 0; w < this->vvdBinCoefs[b].size(); w++)
				{
					double dCoef = this->vvdBinCoefs[b][w];
					if( dCoef )
					{
						double dWavelet = +1.0;
						_ConvetIndexToSub(w, vuSub);

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

				#if	0	// DEL-BY-LEETEN 10/08/2012-BEGIN
					#if	WITH_VECTORS_FOR_COUNTED_COEFS
					for(map<size_t, CWaveletTempCoef>::iterator 
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

					if( !bIsFinalizedWithoutWavelet )
				#endif		// DEL-BY-LEETEN 10/08/2012-END
				for(map<size_t, double>::iterator 
					ipairCoef = this->vmapBinCoefs[b].begin();
					ipairCoef != this->vmapBinCoefs[b].end();
					ipairCoef++)
				{
					double dCoef = ipairCoef->second;
					if( dCoef )
					{
						double dWavelet = 1.0;

						_ConvetIndexToSub(ipairCoef->first, vuSub);

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
				}	// ADD-BY-LEETEN 10/08/2012
			}	
		}

		//! Compute and display statistics for the computed wavelet coefficients.
		virtual	
		void
		_ShowStatistics
		(
		)
		{
			size_t uNrOfNonZeroCoefs = 0;
			for(size_t b = 0; b < UGetNrOfBins(); b++)
			{
				double dEnergy = 0.0;
				for(size_t w = 0; w < this->vvdBinCoefs[b].size(); w++)
				{
					double dCoef = this->vvdBinCoefs[b][w];
					dEnergy += pow(dCoef, 2.0);
					if( dCoef )
						uNrOfNonZeroCoefs++;
				}
				for(map<size_t, double>::iterator
					ipairCoef = vmapBinCoefs[b].begin();
					ipairCoef != vmapBinCoefs[b].end();
					ipairCoef++)
				{
					double dCoef = ipairCoef->second;
					dEnergy += pow(dCoef, 2.0);
					if( dCoef )
						uNrOfNonZeroCoefs++;
				}
				// printf("Energy[%d] = %f\n", b, dEnergy);
			}
			LOG_VAR(uNrOfNonZeroCoefs);

			size_t uNrOfDataItems = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				LOG_VAR(vuDimLengths[d]);
				uNrOfDataItems *= this->vuDataDimLengths[d];
			}
			LOG_VAR(uNrOfDataItems);
			LOG_VAR(UGetNrOfBins());

			double dCR = (double)(uNrOfDataItems * UGetNrOfBins()) / (double)uNrOfNonZeroCoefs;
			LOG_VAR(dCR);

			double dOverhead = (double)uNrOfNonZeroCoefs / (double)uNrOfDataItems;
			LOG_VAR(dOverhead);
		}

		//! Compute statistics of the compressed result.
		virtual
		void
		_SetDimLengths
		(
			const vector<size_t>& vuDimLengths,
			void *_Reserved = NULL
		)
		{
			this->vuDimLengths.clear();
			this->vuDimLevels.clear();

			this->vuDataDimLengths.clear();
			size_t uMaxDimLength = 0;
			for(vector<size_t>::const_iterator 
				ivuDimLength = vuDimLengths.begin();
				ivuDimLength != vuDimLengths.end();
				ivuDimLength++)
			{
				size_t uDimLength = *ivuDimLength;
				this->vuDataDimLengths.push_back(uDimLength);
				uMaxDimLength = max(uMaxDimLength, uDimLength);
			}

			// ADD-BY-LEETEN 10/10/2012-BEGIN
			size_t uMaxDimLevel = (size_t)ceilf(logf((float)uMaxDimLength)/logf(2.0f));
			uMaxDimLength = 1 << uMaxDimLevel;
			// ADD-BY-LEETEN 10/10/2012-END
			for(vector<size_t>::const_iterator 
				ivuDimLength = vuDimLengths.begin();
				ivuDimLength != vuDimLengths.end();
				ivuDimLength++)
			{
				size_t uDimLength = uMaxDimLength;	
				this->vuDimLengths.push_back(uDimLength);
				
				size_t uDimLevel  = (size_t)ceilf(logf((float)uDimLength)/logf(2.0f)) + 1;
				this->vuDimLevels.push_back(uDimLevel);
				vector<size_t> vuSubLevel2Coef;
				vuSubLevel2Coef.resize(uDimLength * uDimLevel);
				for(size_t p = 0, i = 0; i < uDimLength; i++)
				{
					for(size_t 
						l = 0, uW = 1 << uDimLevel, uAccumNrOfWins = 0, uNrOfWins = 1; 
						l < uDimLevel; 
						uAccumNrOfWins += uNrOfWins, uNrOfWins = (!l)?1:(uNrOfWins*2), l++, p++, uW >>= 1)
					{
						size_t uLevelWin = i / uW;
						vuSubLevel2Coef[p] = uAccumNrOfWins + uLevelWin;
					}
				}
				this->vvuSubLevel2Coef.push_back(vuSubLevel2Coef);
			}
			// ADD-BY-LEETEN 10/06/2012-BEGIN
			size_t uThreshold = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				uThreshold *= this->vuDimLengths[d];
			dWaveletThreshold = 1.0 / sqrt((double)uThreshold);
			// ADD-BY-LEETEN 10/06/2012-END
		}
		
		//! Allocate the space to store coefficients for all bins. 
		/*! 
		*/
		virtual	// ADD-BY-LEETEN 09/29/2012
		void 
		_AllocateBins
		(
			size_t uNrOfBins,
			vector<size_t>& vuDimMaxLevels,
			void *_Reserved = NULL
		)
		{
			vdBinWeights.resize(uNrOfBins);	// ADD-BY-LEETEN 10/10/2012

			this->vuDimMaxLevels.clear();
			// MOD-BY-LEETEN 10/19/2012-FROM:	size_t uNrOfCoefs = 1;	// total # coefficients to store
			uNrOfCoefs = 1;
			// MOD-BY-LEETEN 10/19/2012-END

			this->uNrOfUpdatingCoefs = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				// ADD-BY-LEETEN 10/19/2012-BEGIN
				// multiplied by the #coefficients s.t. later the indices can be computed without the extra multiplications
				for(size_t c = 0; c < this->vvuSubLevel2Coef[d].size(); c++)
					this->vvuSubLevel2Coef[d][c] *= uNrOfCoefs;
				// ADD-BY-LEETEN 10/19/2012-END

				size_t uMaxLevel = vuDimMaxLevels[d];
				size_t uCoefLength;	// total number of coefficient for this dimension
				if( !uMaxLevel )
					uMaxLevel = vuDimLevels[d];
				else
					uMaxLevel = min(uMaxLevel, vuDimLevels[d]);
				this->uNrOfUpdatingCoefs *= uMaxLevel;
				uCoefLength = (size_t)1 << (uMaxLevel - 1);
				vuCoefLengths.push_back(uCoefLength);
				this->vuDimMaxLevels.push_back(uMaxLevel);
				uNrOfCoefs *= uCoefLength;		// update the total number of coefficient to store for all dimension
			}

			uNrOfCoefsInFullArray = min(
				(size_t)floor( (double)uSizeOfFullArrays/(double)(uNrOfBins * sizeof(this->vvdBinCoefs[0][0]))), 
				uNrOfCoefs);
			LOG_VAR(uNrOfCoefsInFullArray);

			this->vuCoefDim2Level.resize(uNrOfUpdatingCoefs * UGetNrOfDims());
			for(size_t c = 0; c < uNrOfUpdatingCoefs; c++)
				for(size_t uCoef = c, d = 0; d < UGetNrOfDims(); d++)
				{
					size_t uMaxLevel = this->vuDimMaxLevels[d];
					size_t uLevel = uCoef % uMaxLevel; 
					uCoef /= uMaxLevel; 
					this->vuCoefDim2Level[c * UGetNrOfDims() + d] = uLevel;
 				}

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
			
			// ADD-BY-LEETEN 10/19/2012-BEGIN
			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				size_t uMaxCountPerCoef = 1;	// max # coefficnet
				for(size_t d = 0; d < UGetNrOfDims(); d++, p++)
				{
					size_t uLevel = vuCoefDim2Level[p];
					uMaxCountPerCoef *= (size_t)1 << (vuDimLevels[d] - uLevel);
				}
				vuMaxCounts.push_back(uMaxCountPerCoef);
			}
			// ADD-BY-LEETEN 10/19/2012-END

			#endif	// #if	WITH_VECTORS_FOR_COUNTED_COEFS
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
			// for each bin, apply wavelet transform
			vdSums.clear();

			// ADD-BY-LEETEN 10/06/2012-BEGIN
			size_t uNrOfWavelets = 1;
			for(vector<size_t>::iterator 
				ivuDimMaxLevels = vuDimMaxLevels.begin();
				ivuDimMaxLevels != vuDimMaxLevels.end();
				ivuDimMaxLevels ++)
				uNrOfWavelets *= *ivuDimMaxLevels;

			vector<double> vdWaveletBasis;
			vdWaveletBasis.resize( uNrOfWavelets );
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
				}
			}
			// ADD-BY-LEETEN 10/12/2012-END

			for(size_t b = 0; b < UGetNrOfBins(); b++)
			{
				#if	0	// DEL-BY-LEETEN 10/12/2012-BEGIN
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
					}
				}
				#endif	// DEL-BY-LEETEN 10/12/2012-END

				// now find the combination of the coefficients of all dimensions 
				double dCount = 0.0;	// the combined wavelet basis value. Initially it is one
				for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++, p+= UGetNrOfDims())
				{
					size_t uCoefId = 0;
					#if	0	// MOD-BY-LEETEN 10/19/2012-FROM:
					for(size_t d = 0, uCoefBase = 1;
						d < UGetNrOfDims(); 
						uCoefBase *= vuCoefLengths[d], d++)
					#else		// MOD-BY-LEETEN 10/19/2012-TO:
					for(size_t d = 0; d < UGetNrOfDims(); d++)
					#endif		// MOD-BY-LEETEN 10/19/2012-END
					{
						size_t uLevel = vuCoefDim2Level[p + d];

						// update the index of the current coefficients
						#if	0	// MOD-BY-LEETEN 10/19/2012-FROM:
						size_t uCoef = vvuSubLevel2Coef[d][vuPos[d] * vuDimLevels[d] + uLevel];
						uCoefId += uCoef * uCoefBase;
						#else		// MOD-BY-LEETEN 10/19/2012-TO:
						uCoefId += vvuSubLevel2Coef[d][vuPosLevelProduct[d] + uLevel];
						#endif		// MOD-BY-LEETEN 10/19/2012-END
					}


					double dWaveletCoef = 0.0;
					#if	0	// MOD-BY-LEETEN 10/08/2012-FROM:
					if( uCoefId < this->uNrOfCoefsInFullArray )
						dWaveletCoef = vvdBinCoefs[b][uCoefId];
					else
					{
						map<size_t, double>::iterator ipairCoef = this->vmapBinCoefs[b].find(uCoefId);
						if( this->vmapBinCoefs[b].end() != ipairCoef )
							dWaveletCoef = ipairCoef->second;
					}
					#else		// MOD-BY-LEETEN 10/08/2012-TO:
					dWaveletCoef = DGetBinCoef(b, uCoefId);
					#endif		// MOD-BY-LEETEN 10/08/2012-END

					if( fabs(dWaveletCoef) >= dWaveletThreshold )	// MOD-BY-LEETEN 10/06/2012-FROM:	if( 0.0 != dWaveletCoef )
					{
						double dWavelet = 1.0;
						for(size_t d = 0, uBase = 0;
							d < UGetNrOfDims() && 0.0 != dWavelet; 
							uBase += vuDimMaxLevels[d], d++)
							// combine the wavelet basis value
							dWavelet *= vdWaveletBasis[uBase + vuCoefDim2Level[p + d]];	
							
						// update the corresponding wavelet coeffcients
						if( 0.0 != dWavelet )
							dCount += dWavelet * dWaveletCoef;
					}
				}
				dCount /= dWaveletDenomiator;	
				
				vdSums.push_back(dCount);
			}
		}

		CBase():
			bIsFinalizedWithoutWavelet(false)
		{
		}
	};
}
