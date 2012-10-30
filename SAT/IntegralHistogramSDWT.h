#pragma once

#include <vector>
using namespace std;
#define _USE_MATH_DEFINES
#include <math.h>

// #include "IntegralHistogram.h"
#include "WaveletSAT.h"

template<typename T>
void
_DWT1D
(
	vector<T>& vSrc,
	vector<T>& vDst,
	size_t uLength,
	size_t uLevel,
	void* _Reserved = NULL
)
{
	vector<T>::iterator ivSrcHalf = vSrc.begin();
	for(size_t c = 0; c < uLength/2; c++, ivSrcHalf++)
	{
		T Src1 = vSrc[c * 2];
		T Src2 = vSrc[c * 2 + 1];
		vDst[c] = (Src1 + Src2)/M_SQRT2;
		vDst[c + uLength/2] = (Src1 - Src2)/M_SQRT2;
	}

	if( uLevel > 0 )
	{
		_DWT1D<T>(vDst, vSrc, uLength/2, uLevel - 1);
		std::copy(vSrc.begin(), ivSrcHalf, vDst.begin());
	}
}

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace SAT
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
	template<class T>
	class CIntegralHistogramSDWT:
		public CIntegralHistogram<T>
	{
protected:
		//! Total #coefficients
		size_t uNrOfCoefs;

		//! Product of influencing numbers of coefficients from all dimensions
		/*! 
		#coefficients to update ...
		C = prod l[d]
		*/
		size_t uNrOfUpdatingCoefs;

		//! The maximun number of wavelet levels per dim, which are specified by the user.
		/*!
		m[0], ... m[d], ... m[D - 1]
		*/
		vector<size_t> vuDimMaxLevels;
		
		//! The look up table per dimension to quickly conver dim. subscript to the coefficient indices
		/*!
		For each dim d, 
			the size is l[d] x n[d]: map the spatial location i, i = 0 ... n[d], to its l[d] coefficient locations
		*/
		vector< vector<size_t> > vvuSubLevel2Coef;

		//! #levels per dim
		/*!
		l[0] = log2(n[0]) + 1, ... l[d] = log2(n[d]) + 1, ... l[D - 1] = log2(n[D - 1]) + 1
		*/
		vector<size_t> vuDimLevels;

		//! A lookup table to map the coefficient to its levels per dim.
		/*! size: D x C: 
		Map each coefficients c, c = 0 ... C - 1 = prod l[d] - 1, to the the corresponding D levels
		*/
		vector<size_t> vuCoefDim2Level;

		//! The denomator for Wavelet basis
		/*!
		W = sqrt(prod n[0], ... n[d], ... n[D])
		*/
		double dWaveletDenomiator;

		//! The mininal value of the projected wavelet coefficients
		double dWaveletThreshold;	

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
			dCoef = this->vvdBinSATs[uBin][uCoefId];
			return dCoef;
		}
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

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

		//! Compute statistics of the compressed result.
		virtual
		void
		_SetDimLengths
		(
			const vector<size_t>& vuDimLengths,
			void *_Reserved = NULL
		)
		{
			CIntegralHistogram<T>::_SetDimLengths(vuDimLengths);

			this->vuDimLengths.clear();
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
			size_t uThreshold = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				uThreshold *= this->vuDimLengths[d];
			dWaveletThreshold = 1.0 / sqrt((double)uThreshold);
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
			CIntegralHistogram<T>::_AllocateBins(uNrOfBins, vuDimMaxLevels);

			this->vuDimMaxLevels.clear();
			uNrOfCoefs = 1;

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
				this->vuDimMaxLevels.push_back(uMaxLevel);
				uNrOfCoefs *= uCoefLength;		// update the total number of coefficient to store for all dimension
			}

			this->vuCoefDim2Level.resize(uNrOfUpdatingCoefs * UGetNrOfDims());
			for(size_t c = 0; c < uNrOfUpdatingCoefs; c++)
				for(size_t uCoef = c, d = 0; d < UGetNrOfDims(); d++)
				{
					size_t uMaxLevel = this->vuDimMaxLevels[d];
					size_t uLevel = uCoef % uMaxLevel; 
					uCoef /= uMaxLevel; 
					this->vuCoefDim2Level[c * UGetNrOfDims() + d] = uLevel;
 				}
		}

		//! Finalize the computation of SAT
		virtual 
		void 
		_Finalize
		(
			void *_Reserved = NULL
		)
		{
			CIntegralHistogram<T>::_Finalize();
			
			int iWaveletDenomiator = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				iWaveletDenomiator *= 1 << (vuDimLevels[d] - 1);
			dWaveletDenomiator = sqrt((double)iWaveletDenomiator);

			vector<size_t> vuScanLineBase;
			vuScanLineBase.resize(this->UGetNrOfDims());

			for(size_t uBase = 1, d = 0; 
				d < this->UGetNrOfDims(); 
				uBase *= this->vuDimLengths[d], d++)
			{
				size_t uNrOfLevels;
				uNrOfLevels = (size_t)ceil(log((double)this->vuDimLengths[d])/log(2.0));

				size_t uNrOfScanLines = this->uDataLength / this->vuDimLengths[d];

				vector<size_t> vuScanLineIndices;
				vuScanLineIndices.resize(this->vuDimLengths[d]);

				vector<T> vSrc;
				vSrc.resize(this->vuDimLengths[d]);

				vector<T> vDst;
				vDst.resize(this->vuDimLengths[d]);

				for(size_t i = 0; i < uNrOfScanLines; i++)
				{
					for(size_t uIndex = i, d2 = 0; d2 < this->UGetNrOfDims(); d2++)
					{
						if( d2 == d )
							vuScanLineBase[d2] = 0;
						else
						{
							vuScanLineBase[d2] = (uIndex % this->vuDimLengths[d2]);
							uIndex /= this->vuDimLengths[d2];
						}
					}

					// store the element indices along the current scanline so it can be reused for all bins
					vuScanLineIndices[0] = this->UConvertSubToIndex(vuScanLineBase);
					for(size_t j = 1; j < this->vuDimLengths[d]; j++)
					{
						vuScanLineIndices[j] = vuScanLineIndices[j - 1] + uBase;
					}

					for(size_t b = 0; b < UGetNrOfBins(); b++)
					{
						for(size_t j = 0; j < this->vuDimLengths[d]; j++)
							vSrc[j] = this->vvdBinSATs[b][vuScanLineIndices[j]];
						_DWT1D(vSrc, vDst, this->vuDimLengths[d], uNrOfLevels - 1);
						for(size_t j = 0; j < this->vuDimLengths[d]; j++)
							this->vvdBinSATs[b][vuScanLineIndices[j]] = vDst[j];
					}
				}
			}
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

			size_t uNrOfWavelets = 0;
			for(vector<size_t>::iterator 
				ivuDimMaxLevels = vuDimMaxLevels.begin();
				ivuDimMaxLevels != vuDimMaxLevels.end();
				ivuDimMaxLevels ++)
				uNrOfWavelets += *ivuDimMaxLevels;

			vector<double> vdWaveletBasis;
			vdWaveletBasis.resize( uNrOfWavelets );

			vector<size_t> vuPosLevelProduct;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				vuPosLevelProduct.push_back(vuPos[d] * vuDimLevels[d]);

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
		
			vdSums.resize(UGetNrOfBins());

			// now find the combination of the coefficients of all dimensions 
			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				size_t uCoefId = 0;
				double dWavelet = 1.0;
				for(size_t d = 0, uBase = 0;
					d < UGetNrOfDims(); 
					uBase += vuDimMaxLevels[d], d++, p++)
				{
					// update the index of the current coefficients
					uCoefId += vvuSubLevel2Coef[d][vuPosLevelProduct[d] + vuCoefDim2Level[p]];

					// combine the wavelet basis value
					dWavelet *= vdWaveletBasis[uBase + vuCoefDim2Level[p]];	
				}
				for(size_t b = 0; b < UGetNrOfBins(); b++)
				{
					double dWaveletCoef = DGetBinCoef(b, uCoefId);

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

		CIntegralHistogramSDWT()
		{
		}
	};
}
