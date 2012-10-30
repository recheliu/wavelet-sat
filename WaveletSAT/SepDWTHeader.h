#pragma once

//! Decide whether a table of the product of wavelet bais is precomputed.
/*!
This table is used in the decompression stage.
*/
#define	WITH_PRECOMPUTED_WAVELET_BASIS	0	

//! Decide whether a table of the product of wavelet sum is precomputed.
/*!
This table is used in the compression stage.
*/
#define	WITH_PRECOMPUTED_WAVELET_SUMS	0	

//! Decide whether a table that immediately map each updating coefficients and its dimension to the corresponding 1D index in the wavelet table per dimension.
#define WITH_COEF_DIM_2_WAVELET		1
// ADD-BY-LEETEN 10/21/2012-END

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

#include "HeaderBase.h"	
#include "DWT.h"
#include "liblog.h"	

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
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
	class CSepDWTHeader:
		public CHeaderBase	// ADD-BY-LEETEN 09/29/2012
	{
protected:	
		// ADD-BY-LEETEN 10/19/2012-BEGIN
		//! Total #coefficients
		size_t uNrOfCoefs;
		// ADD-BY-LEETEN 10/19/2012-END

		//! #Coefs per dim.
		/*!
		c[0], ..., c[d], ... c[D - 1]
		*/
		vector<size_t> vuCoefLengths;

		//! The weight sum for each bin
		vector<double> vdBinWeights;

		//! The maximun number of wavelet levels per dim, which are specified by the user.
		/*!
		m[0], ... m[d], ... m[D - 1]
		*/
		vector<size_t> vuDimMaxLevels;
		
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

		// ADD-BY-LEETEN 10/18/2012-BEGIN
		#if	WITH_PRECOMPUTED_WAVELET_BASIS
		//! The volume to store the wavelet basis for all updating coefficients.
		vector<double>	vdWaveletBasisPerUpdatingCoefs;
		#endif	// #if	WITH_PRECOMPUTED_WAVELET_BASIS
		// ADD-BY-LEETEN 10/18/2012-END

		// ADD-BY-LEETEN 10/21/2012-BEGIN
		#if	WITH_PRECOMPUTED_WAVELET_SUMS
		vector<double>	vdWaveletSigns;
		vector<double>	vdWaveletSums;
		#endif	// if	WITH_PRECOMPUTED_WAVELET_SUMS
		// ADD-BY-LEETEN 10/21/2012-END

		//! A lookup table to map the coefficient to its levels per dim.
		/*! size: D x C: 
		Map each coefficients c, c = 0 ... C - 1 = prod l[d] - 1, to the the corresponding D levels
		*/
		vector<size_t> vuCoefDim2Level;

		// ADD-BY-LEETEN 10/21/2012-BEGIN
		#if	WITH_COEF_DIM_2_WAVELET
		//! A lookup table to map the coefficient to its 1D index in the L x D wavelets
		/*! size: D x C: 
		Map each coefficients c, c = 0 ... C - 1 = prod l[d] - 1, to the the corresponding index 
		*/
		vector<size_t> vuCoefDim2Wavelet;
		#endif	// #if	WITH_COEF_DIM_2_WAVELET
		// ADD-BY-LEETEN 10/21/2012-END

		//! The look up table per dimension to quickly conver dim. subscript to the coefficient indices
		/*!
		For each dim d, 
			the size is l[d] x n[d]: map the spatial location i, i = 0 ... n[d], to its l[d] coefficient locations
		*/
		vector< vector<size_t> > vvuSubLevel2Coef;

		//! Max. count per coefficients
		vector<size_t> vuMaxCounts;	// ADD-BY-LEETEN 10/19/2012

		//! The denomator for Wavelet basis
		/*!
		W = sqrt(prod n[0], ... n[d], ... n[D])
		*/
		double dWaveletDenomiator;

		//! The mininal value of the projected wavelet coefficients
		double dWaveletThreshold;	
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0100,
			PARAMETER_END
		};

		//! Retunr the smallest value. 
		/*!
		This can be used to filter too small value caused by numerical error.
		*/
		virtual	
		double 
		DGetThreshold
		(
			void *_Reserved = NULL
		)
		{
			return dWaveletThreshold;
		}

		virtual	
		void 
		_SetMaxLevels
		(
			vector<size_t>& vuDimMaxLevels,
			void *_Reserved = NULL
		)
		{
			this->vuDimMaxLevels.clear();
			uNrOfCoefs = 1;
			this->uNrOfUpdatingCoefs = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				#if	0	// DEL-BY-LEETEN 10/29/2012-BEGIN
				// multiplied by the #coefficients s.t. later the indices can be computed without the extra multiplications
				for(size_t c = 0; c < this->vvuSubLevel2Coef[d].size(); c++)
					this->vvuSubLevel2Coef[d][c] *= uNrOfCoefs;
				#endif		// DEL-BY-LEETEN 10/29/2012-END

				size_t uMaxLevel = vuDimMaxLevels[d];
				size_t uCoefLength;	// total number of coefficient for this dimension
				if( !uMaxLevel )
					uMaxLevel = vuDimLevels[d];
				else
					uMaxLevel = min(uMaxLevel, vuDimLevels[d]);
				this->uNrOfUpdatingCoefs *= uMaxLevel;
				uCoefLength = (size_t)1 << (uMaxLevel - 1);
				vuCoefLengths[d] = uCoefLength;	// MOD-BY-LEETEN 10/29/2012-FROM:	vuCoefLengths.push_back(uCoefLength);
				this->vuDimMaxLevels.push_back(uMaxLevel);
				uNrOfCoefs *= uCoefLength;		// update the total number of coefficient to store for all dimension
			}

			this->vuCoefDim2Level.resize(uNrOfUpdatingCoefs * UGetNrOfDims());
			#if	WITH_COEF_DIM_2_WAVELET			
			this->vuCoefDim2Wavelet.resize(uNrOfUpdatingCoefs * UGetNrOfDims());
			#endif	// #if	WITH_COEF_DIM_2_WAVELET		
			for(size_t c = 0; c < uNrOfUpdatingCoefs; c++)
				#if	!WITH_COEF_DIM_2_WAVELET		// ADD-BY-LEETEN 10/21/2012
				for(size_t uCoef = c, d = 0; d < UGetNrOfDims(); d++)
				// ADD-BY-LEETEN 10/21/2012-BEGIN
				#else	// #if	!WITH_COEF_DIM_2_WAVELET	
				for(size_t	uCoef = c, uBase = 0,	d = 0; 
						d < UGetNrOfDims(); 
						uBase += this->vuDimMaxLevels[d], d++)
				#endif	// #if	!WITH_COEF_DIM_2_WAVELET	
				// ADD-BY-LEETEN 10/21/2012-END
				{
					size_t uMaxLevel = this->vuDimMaxLevels[d];
					size_t uLevel = uCoef % uMaxLevel; 
					uCoef /= uMaxLevel; 
					this->vuCoefDim2Level[c * UGetNrOfDims() + d] = uLevel;
					#if	WITH_COEF_DIM_2_WAVELET		
					this->vuCoefDim2Wavelet[c * UGetNrOfDims() + d] = uBase + uLevel;
					#endif	// #if	WITH_COEF_DIM_2_WAVELET	
 				}

			// ADD-BY-LEETEN 10/18/2012-BEGIN
			#if	WITH_PRECOMPUTED_WAVELET_BASIS
			vector<double> vdWaveletBasis;
			vdWaveletBasis.resize( this->UGetNrOfDims() * this->vuDimMaxLevels[0] );

			// for each dimension d, based on the posistion, store the corresponding l[d] wavelet basis value
			for(size_t p = 0,	d = 0; d < this->UGetNrOfDims();	d++)
				for(size_t	l = 0; l < this->vuDimMaxLevels[d];	l++, p++)
					vdWaveletBasis[p] = ( l >= 2 )?sqrt( (double)(1 << (l - 1)) ):1.0;

			vdWaveletBasisPerUpdatingCoefs.resize(uNrOfUpdatingCoefs);
			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				double dWavelet = 1.0;
				for(size_t d = 0, uBase = 0;
					d < UGetNrOfDims(); 
					uBase += this->vuDimMaxLevels[d], d++, p++)
					// combine the wavelet basis value
					dWavelet *= vdWaveletBasis[uBase + vuCoefDim2Level[p]];	
				vdWaveletBasisPerUpdatingCoefs[c] = dWavelet;
			}
			#endif	// #if	WITH_PRECOMPUTED_WAVELET_BASIS

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

			#if	WITH_PRECOMPUTED_WAVELET_SUMS
			vdWaveletSigns.resize(uNrOfUpdatingCoefs);
			for(size_t p = 0,	c = 0; c < uNrOfUpdatingCoefs;	c++)
			{
				int iSign = +1;
				for(size_t	d = 0; d < UGetNrOfDims();	d++, p++)
					if( this->vuCoefDim2Level[p] )
						iSign = -iSign;
				vdWaveletSigns[c] = (double)iSign;
			}

			size_t uWaveletSumSize = 1;
			for(size_t d = 0; d < this->UGetNrOfDims(); d++)
				uWaveletSumSize *= this->vuDimLengths[d] + 1;	
			vdWaveletSums.resize(uWaveletSumSize);

			for(size_t c = 0; c < uWaveletSumSize; c++)
			{
				long lProd = 1;
				for(size_t uIndex = c, d = 0; d < UGetNrOfDims(); d++)
				{
					lProd *= uIndex % (this->vuDimLengths[d] + 1);
					uIndex /= this->vuDimLengths[d] + 1;
				}
				vdWaveletSums[c] = (double)lProd;
			}
			#endif	// if	WITH_PRECOMPUTED_WAVELET_SUMS
		}

		//! Compute statistics of the compressed result.
		virtual
		void
		_Set
		(
			const vector<size_t>& vuDimLengths,
			const size_t uNrOfBins,
			void *_Reserved = NULL
		)
		{
			CHeaderBase::_Set(vuDimLengths, uNrOfBins);

			vdBinWeights.resize(uNrOfBins);

			this->vuCoefLengths.clear();
			this->vuDimLevels.clear();

			size_t uMaxDimLength = 0;
			for(vector<size_t>::const_iterator 
				ivuDimLength = vuDimLengths.begin();
				ivuDimLength != vuDimLengths.end();
				ivuDimLength++)
			{
				size_t uDimLength = *ivuDimLength;
				uMaxDimLength = max(uMaxDimLength, uDimLength);
			}

			size_t uMaxDimLevel = (size_t)ceilf(logf((float)uMaxDimLength)/logf(2.0f));
			uMaxDimLength = (size_t)1 << uMaxDimLevel;
			for(vector<size_t>::const_iterator 
				ivuDimLength = vuDimLengths.begin();
				ivuDimLength != vuDimLengths.end();
				ivuDimLength++)
			{
				size_t uCoefLength = uMaxDimLength;	
				this->vuCoefLengths.push_back(uCoefLength);
				size_t uDimLevel  = (size_t)ceilf(logf((float)uCoefLength)/logf(2.0f)) + 1;
				this->vuDimLevels.push_back(uDimLevel);
				vector<size_t> vuSubLevel2Coef;
				vuSubLevel2Coef.resize(uCoefLength * uDimLevel);
				for(size_t p = 0, i = 0; i < uCoefLength; i++)
					for(size_t 
						l = 0, uW = 1 << uDimLevel, uAccumNrOfWins = 0, uNrOfWins = 1; 
						l < uDimLevel; 
						uAccumNrOfWins += uNrOfWins, uNrOfWins = (!l)?1:(uNrOfWins*2), l++, p++, uW >>= 1)
					{
						size_t uLevelWin = i / uW;
						vuSubLevel2Coef[p] = uAccumNrOfWins + uLevelWin;
					}

				this->vvuSubLevel2Coef.push_back(vuSubLevel2Coef);
			}

			int iWaveletDenomiator = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				iWaveletDenomiator *= 1 << (vuDimLevels[d] - 1);
			dWaveletDenomiator = sqrt((double)iWaveletDenomiator);

			size_t uThreshold = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				uThreshold *= this->vuCoefLengths[d];
			dWaveletThreshold = 1.0 / sqrt((double)uThreshold);

			_SetMaxLevels(this->vuDimLevels);
		}
		
		CSepDWTHeader()
		{
		}
	};
}
