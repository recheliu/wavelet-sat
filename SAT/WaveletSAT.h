#pragma once

#include <vector>
using namespace std;
#include <math.h>

#include "liblog.h"	// ADD-BY-LEETEN	09/09/2012

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

	template<class T>
	class CBase
	{
		//! The maximun number of wavelet levels per dim. 
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
		
		//!
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
		size_t UGetNrOfDims() const {	return vuDimLengths.size();	};
		
		//! #Dimensions
		/*!
		B
		*/
		size_t UGetNrOfBins() const {	return vvdBinCoefs.size();	};

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
		// MOD-BY-LEETEN 09/09/2012-FROM: vector<vector<size_t>> vvuSubLevel2Coef;
		vector< vector<size_t> > vvuSubLevel2Coef;
		// MOD-BY-LEETEN 09/09/2012-END

		// ADD-BY-LEETEN 09/07/2012-BEGIN
		//! The denomator for Wavelet basis
		/*!
		W = sqrt(prod n[0], ... n[d], ... n[D])
		*/
		double dWaveletDenomiator;
		// ADD-BY-LEETEN 09/07/2012-END
protected:		
		void 
		_UpdateBin
		(
			const vector<size_t>& vuPos, 
			const T& value,
			size_t uBin, 
			double dWeight,	// ADD-BY-LEETEN 09/07/2012
			void *_Reserved = NULL
		)
		{
			vector<double> vdWavelets;
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
					double dWavelet;
					
					// Given the currenet level l and subscript uPos, compute the sum of the portion in wavelet after uPos
					size_t uPosInWavelet = uPos % w;
					if( 0 == l )
						dWavelet = (double)(w / 2 - uPosInWavelet);
					else
					{
						if( uPosInWavelet < w / 2)
							dWavelet = (double)uPosInWavelet;
						else
							dWavelet = (double)(w - uPosInWavelet);
						dWavelet *= -sqrt((double)(1 << (l - 1) ));	
					}
					// DEL-BY-LEETEN 09/07/2012:	dWavelet /= sqrt((double)(uMaxWin/2));
					
					vdWavelets.push_back( dWavelet );
				}
			}
			
			// compute the weight from this value
			#if	0	// DEL-BY-LEETEN 09/07/2012-BEGIN
				double dW;
				_MapValueToBinWeight(vuPos, value, uBin, dW);
			#endif		// DEL-BY-LEETEN 09/07/2012-END

			// now find the combination of the coefficients of all dimensions
			for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				// MOD-BY-LEETEN 09/07/2012-FROM:	double dWavelet = dW;
				double dWavelet = dWeight;
				// MOD-BY-LEETEN 09/07/2012-END

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
					
					dWavelet *= vdWavelets[uBase + uLevel];
					
					size_t uCoef = vvuSubLevel2Coef[d][vuPos[d] * vuDimLevels[d] + uLevel];
					uCoefId += uCoef * uCoefBase;
				}
				// update the corresponding wavelet coeffcients
				vvdBinCoefs[uBin][uCoefId] += dWavelet;
			}
		}

		//! Given the subscript of a point, update the corresponding bin SAT
		/*!
		Note: As this method is not public, it is not intended to be directly used by the application. 
		Instead, the sub class should define additional methods to call this method for each all grid points.
		This design is for the case that not all data's area available
		*/
		void 
		_Update
		(
			const vector<size_t>& vuPos,
			const T& value,
			void *_Reserved = NULL
		)
		{
			#if	0	// MOD-BY-LEETEN 09/07/2012-FROM:
				vector<size_t> vuBins;
				_MapValueToBins(vuPos, value, vuBins);
				for(vector<size_t>::iterator ivuBin  = vuBins.begin(); ivuBin != vuBins.end(); ivuBin++)
					_UpdateBin(vuPos, value, *ivuBin);
			#else		// MOD-BY-LEETEN 09/07/2012-TO:
			// MOD-BY-LEETEN 09/09/2012-FROM: vector<pair<size_t, double>> vpBins;
			vector< pair<size_t, double> > vpBins;
			// MOD-BY-LEETEN 09/09/2012-END

			_MapValueToBins(vuPos, value, vpBins);

#if 0 // MOD-BY-LEETEN 09/09/2012-FROM:
			for(vector<pair<size_t, double>>::iterator ivpBin  = vpBins.begin(); ivpBin != vpBins.end(); ivpBin++)
#else // MOD-BY-LEETEN 09/09/2012-TO:
			for(vector< pair<size_t, double> >::iterator ivpBin  = vpBins.begin(); ivpBin != vpBins.end(); ivpBin++)
#endif // MOD-BY-LEETEN 09/09/2012-END
				_UpdateBin(vuPos, value, ivpBin->first, ivpBin->second);
			#endif		// MOD-BY-LEETEN 09/07/2012-END
			
		}

		////////////////////////////////////////////////////////////////////
		/*
		The interface that should be overloaed by the sub class
		*/

		#if	0	// DEL-BY-LEETEN 09/07/2012-BEGIN
			//! Map the given value to the weight contributed to the given bin in the histogram
			virtual
			void
			_MapValueToBinWeight
			(
				const vector<size_t>& vuPos,
				const T& value, 
				size_t uBin,
				double& dW,
				void *_Reserved = NULL
			) = 0;
		#endif		// DEL-BY-LEETEN 09/07/2012-END
		
		//! This method should be overloaded to return the indices of bins that will be changed
		virtual 
		void 
		_MapValueToBins
		(
			const vector<size_t>& vuPos,
			const T& value, 
			// MOD-BY-LEETEN 09/07/2012-FROM:			vector<size_t>& vuBins,
			// MOD-BY-LEETEN 09/09/2012-FROM: vector<pair<size_t, double>>& vuBins,
			vector< pair<size_t, double> >& vuBins,
			// MOD-BY-LEETEN 09/09/2012-END
			// MOD-BY-LEETEN 09/07/2012-END
			void *_Reserved = NULL
		) = 0;
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/
		
		// ADD-BY-LEETEN 09/07/2012-BEGIN
		//! Finalize the computation of SAT
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

			for(size_t b = 0; b < UGetNrOfBins(); b++)
				for(size_t w = 0; w < this->vvdBinCoefs[b].size(); w++)
				{
					double dCoef = this->vvdBinCoefs[b][w];
					if( dCoef )
						this->vvdBinCoefs[b][w] /= dWaveletDenomiator;
				}
		}
		// ADD-BY-LEETEN 09/07/2012-END

		#if	0	// MOD-BY-LEETEN	09/09/2012-FROM:
		//! Get the energy of all bin SATs.
		void
		_GetEnergy
		(
			vector<double>& vdBinEnergies,
			void *_Reserved = NULL
		)
		{
			vdBinEnergies.clear();
			for(size_t b = 0; b < UGetNrOfBins(); b++)
			{
				double dEnergy = 0.0;
				for(size_t w = 0; w < this->vvdBinCoefs[b].size(); w++)
				{
					double dCoef = this->vvdBinCoefs[b][w];
					dEnergy += pow(dCoef, 2.0);
				}
				vdBinEnergies.push_back(dEnergy);
			}
		}
		#else	// MOD-BY-LEETEN	09/09/2012-TO:
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
				// printf("Energy[%d] = %f\n", b, dEnergy);
			}
			LOG_VAR(uNrOfNonZeroCoefs);

			size_t uNrOfDataItems = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				LOG_VAR(vuDimLengths[d]);
				uNrOfDataItems *= this->vuDimLengths[d];
			}
			LOG_VAR(uNrOfDataItems);
			LOG_VAR(UGetNrOfBins());

			double dCR = (double)(uNrOfDataItems * UGetNrOfBins()) / (double)uNrOfNonZeroCoefs;
			LOG_VAR(dCR);

			double dOverhead = (double)uNrOfNonZeroCoefs / (double)uNrOfDataItems;
			LOG_VAR(dOverhead);
		}
		#endif	// MOD-BY-LEETEN	09/09/2012-END

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
			for(vector<size_t>::const_iterator 
				ivuDimLength = vuDimLengths.begin();
				ivuDimLength != vuDimLengths.end();
				ivuDimLength++)
			{
				size_t uDimLength = *ivuDimLength;
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
		}
		
		//! Allocate the space to store coefficients for all bins. 
		/*! 
		*/
		void 
		_AllocateBins
		(
			size_t uNrOfBins,
			vector<size_t>& vuDimMaxLevels,
			void *_Reserved = NULL
		)
		{
			this->vuDimMaxLevels.clear();
			size_t uNrOfCoefs = 1;	// total # coefficients to store
			this->uNrOfUpdatingCoefs = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				size_t uMaxLevel = vuDimMaxLevels[d];
				size_t uCoefLength;	// total number of coefficient for this dimension
				if( !uMaxLevel )
					uMaxLevel = vuDimLevels[d];
				else
					uMaxLevel = min(uMaxLevel, vuDimLevels[d]);
				this->uNrOfUpdatingCoefs *= uMaxLevel;
				// MOD-BY-LEETEN 09/09/2012-FROM:				uCoefLength = 1 << (uMaxLevel - 1);
				uCoefLength = (size_t)1 << (uMaxLevel - 1);
				// MOD-BY-LEETEN 09/09/2012-END
				vuCoefLengths.push_back(uCoefLength);
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

			vvdBinCoefs.resize(uNrOfBins);
			for(size_t b = 0; b < uNrOfBins; b++)
				vvdBinCoefs[b].resize(uNrOfCoefs);
		}
		
		//! Return the sum of all bins at the given position
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
			for(size_t b = 0; b < UGetNrOfBins(); b++)
			{
				// for each dimenion d, based on the posistion, store the corresponding l[d] wavelet basis value

				vector<double> vdWaveletBasis;
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
						// DEL-BY-LEETEN 09/07/2012:	dWaveletBasis /= sqrt( (double)(uMaxWin/2) );
					
						vdWaveletBasis.push_back( dWaveletBasis );
					}
				}
				// now find the combination of the coefficients of all dimensions 
				double dCount = 0.0;	// the combined wavelet basis value. Initially it is one
				for(size_t p = 0, c = 0; c < uNrOfUpdatingCoefs; c++, p+= UGetNrOfDims())
				{
					size_t uCoefId = 0;
					for(size_t d = 0, uCoefBase = 1;
						d < UGetNrOfDims(); 
						uCoefBase *= vuCoefLengths[d], d++)
					{
						size_t uLevel = vuCoefDim2Level[p + d];

						// update the index of the current coefficients
						size_t uCoef = vvuSubLevel2Coef[d][vuPos[d] * vuDimLevels[d] + uLevel];
						uCoefId += uCoef * uCoefBase;
					}

					if( 0.0 != vvdBinCoefs[b][uCoefId] )
					{
						double dWavelet = 1.0;
						for(size_t d = 0, uBase = 0;
							d < UGetNrOfDims() && 0.0 != dWavelet; 
							uBase += vuDimMaxLevels[d], d++)
						{
							size_t uLevel = vuCoefDim2Level[p + d];
							
							// combine the wavelet basis value
							dWavelet *= vdWaveletBasis[uBase + uLevel];	
						}
						
						// update the corresponding wavelet coeffcients
						if( 0.0 != dWavelet )
							dCount += dWavelet * vvdBinCoefs[b][uCoefId];
					}
				}
				dCount /= dWaveletDenomiator;	// ADD-BY-LEETEN 09/07/2012
				vdSums.push_back(dCount);
			}
		}
	};
}
