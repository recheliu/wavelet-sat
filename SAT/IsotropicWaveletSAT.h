#pragma once

// ADD-BY-LEETEN 09/12/2012-BEGIN

#define WITH_VECTORS_FOR_COUNTED_COEFS	1	// ADD-BY-LEETEN 09/14/2012

#include <vector>
using namespace std;
#include <math.h>

#include "WaveletSAT.h"	// ADD-BY-LEETEN 09/29/2012

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	template<class T>
	class CIsotropicWaveletSAT
		:public CBase<T>	// ADD-BY-LEETEN 09/29/2012
	{
	// ADD-BY-LEETEN 10/01/2012-BEGIN
	protected:
		// DEL-BY-LEETEN 10/08/2012:	vector< vector<double> > vvdBinIsotropicCoefs;

	// ADD-BY-LEETEN 10/05/2012-BEGIN
	//! #Wavelet per level. 
	/*!
	W. The value should be 2^D - 1.
	*/
	size_t uNrOfWaveletsPerLevel;

	//! A table that maps the W wavelet to the corresponding binary offset
	/*
	For each of the W wavelets, this table maps them to the corresponding binary offset.
	For instance, for 2D case, which has 3 wavelets, this table has 3 x 2 entries:
		0 1
		1 0
		1 1 

	For 3D, the table has 7 x 3 entries:
		0 0 1
		0 1 0
		0 1 1
		1 0 0
		1 0 1
		1 1 0
		1 1 1
	*/
	vector< size_t > vuWaveletToOffset;
	// ADD-BY-LEETEN 10/05/2012-END

	// ADD-BY-LEETEN 10/01/2012-END
public:
		virtual	// ADD-BY-LEETEN 09/29/2012
		void
		_SetLong(
			int eName,
			long lValue,
			void* _Reserved = NULL
		)
		{
			switch(eName)
			{
			default:
				;
			}
			CBase<T>::_SetLong(eName, lValue);
		}
		
		#if	0	// DEL-BY-LEETEN 10/08/2012-BEGIN
			// ADD-BY-LEETEN 10/01/2012-BEGIN
			virtual	
			void
			_ShowStatistics
			(
			)
			{
			  CBase<T>::_ShowStatistics();
				size_t uNrOfNonZeroCoefs = 0;
				for(size_t b = 0; b < this->UGetNrOfBins(); b++)
				{
					double dEnergy = 0.0;
					for(size_t w = 0; w < this->vvdBinIsotropicCoefs[b].size(); w++)
					{
						double dCoef = this->vvdBinIsotropicCoefs[b][w];
						dEnergy += pow(dCoef, 2.0);
						if( dCoef )
							uNrOfNonZeroCoefs++;
					}
					/*
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
					*/
				}
				LOG_VAR(uNrOfNonZeroCoefs);

				size_t uNrOfDataItems = 1;
				for(size_t d = 0; d < this->UGetNrOfDims(); d++)
				{
					LOG_VAR(this->vuDimLengths[d]);
					uNrOfDataItems *= this->vuDataDimLengths[d];
				}
				LOG_VAR(uNrOfDataItems);
				LOG_VAR(this->UGetNrOfBins());

				double dCR = (double)(uNrOfDataItems * this->UGetNrOfBins()) / (double)uNrOfNonZeroCoefs;
				LOG_VAR(dCR);

				double dOverhead = (double)uNrOfNonZeroCoefs / (double)uNrOfDataItems;
				LOG_VAR(dOverhead);
			}
			// ADD-BY-LEETEN 10/01/2012-END
		#endif	// DEL-BY-LEETEN 10/08/2012-END

		// ADD-BY-LEETEN 10/05/2012-BEGIN
		virtual
		void
		_SetDimLengths
		(
			const vector<size_t>& vuDimLengths,
			void *_Reserved = NULL
		)
		{
			CBase<T>::_SetDimLengths(vuDimLengths);

			// Decide W = 2^D - 1
			uNrOfWaveletsPerLevel = ((size_t)1 << this->UGetNrOfDims()) - 1;

			// build the table that maps each wavelet to its binary representation
			vuWaveletToOffset.clear();
			for(size_t w = 0; w < uNrOfWaveletsPerLevel; w++)
				for(size_t i = w + 1, d = 0; d < this->UGetNrOfDims(); i /= 2, d++)
					vuWaveletToOffset.push_back( i % 2 );
		}
		// ADD-BY-LEETEN 10/05/2012-END

		//! Return the sum of all bins at the given position
		virtual	// ADD-BY-LEETEN 09/29/2012
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
			vector<size_t> vuWaveletPos;
			vuWaveletPos.resize(vuPos.size());

			vector<size_t> vuPosInWavelet;
			vuPosInWavelet.resize(vuPos.size());

			vector< size_t > vuCoefPos;
			vuCoefPos.resize(vuPos.size());
			// ADD-BY-LEETEN 10/06/2012-END

			for(size_t b = 0; b < this->UGetNrOfBins(); b++) // MOD-BY-LEETEN 10/05/2012: for(size_t b = 0; b < UGetNrOfBins(); b++)
			{
				double dCount = 0.0;

				// get the scale coefficients
				// Here we assume Harr transform
				dCount += this->DGetBinCoef(b, 0);	// MOD-BY-LEETEN 10/08/2012-FROM:	dCount += this->vvdBinIsotropicCoefs[b][0];

				// now get the wavelet coefficients
				for(size_t 
					l = 0, uWaveletSize = (size_t)1 << (this->vuDimLevels[0] - 1);
					l < (this->vuDimLevels[0] - 1); 
					l++, uWaveletSize >>= 1)
				{
					double dWaveletWeight = pow( (double)(1 << l), (double)this->UGetNrOfDims()/2.0);

					// get the wavelet position
					for(size_t d = 0; d < vuPos.size(); d++)
					{
						vuWaveletPos[d] = vuPos[d] / uWaveletSize;
						vuPosInWavelet[d] = vuPos[d] % uWaveletSize;
					}

					double dWaveletMaxAbsCoef = 0.0;	// ADD-BY-LEETEN 10/06/2012

					for(size_t p = 0, w = 0; w < uNrOfWaveletsPerLevel; w++, p += this->UGetNrOfDims())
					{
						// locate the wavelet coefficients
						for(size_t d = 0; d < this->UGetNrOfDims(); d++)
							vuCoefPos[d] = vuWaveletToOffset[p + d] * ((size_t) 1 << l) + vuWaveletPos[d];

						size_t uIndex = this->UConvetSubToIndex(vuCoefPos);
						double dWaveletCoef = this->DGetBinCoef(b, uIndex);	// MOD-BY-LEETEN 10/08/2012-FROM:	double dWaveletCoef = this->vvdBinIsotropicCoefs[b][uIndex];

						// now find the basis
						if( fabs(dWaveletCoef) >= this->dWaveletThreshold )
						{
							// here assume that it is Harr transform
							// decide the sign
							double dWavelet = dWaveletWeight;
							for(size_t d = 0; d < this->UGetNrOfDims(); d++)
								if( vuWaveletToOffset[p + d] && vuPosInWavelet[d] >= uWaveletSize / 2)
									dWavelet = -dWavelet;

							dCount += dWaveletCoef * dWavelet;

							dWaveletMaxAbsCoef = max(dWaveletMaxAbsCoef, (double)fabs(dWaveletCoef ));	// ADD-BY-LEETEN 10/06/2012
						}
					}
					// ADD-BY-LEETEN 10/06/2012-BEGIN
					if( !dWaveletMaxAbsCoef )
						break;
					// ADD-BY-LEETEN 10/06/2012-END
				}
				if( dCount )
					dCount /= this->dWaveletDenomiator;

				vdSums.push_back(dCount);
			} // for b
		}
		
		//! Finalize the computation of SAT
		virtual	
		void 
		_Finalize
		(
			void *_Reserved = NULL
		)
		{
			// First, finish the computation of the (separable) DWT
			// Do not apply wavelet coefficients
			this->_SetBoolean(CBase<T>::FINALIZED_WITHOUT_WAVELET, true);
			CBase<T>::_Finalize
			(
			);

			// ADD-BY-LEETEN 10/06/2012-BEGIN
			vector<size_t> vuBase;
			vuBase.resize(this->UGetNrOfDims());

			vector<size_t> vuSrc;
			vuSrc.resize(this->UGetNrOfDims()); 

			vector<double> vdSrc;
			vdSrc.resize(this->vuDimLengths[0]);

			vector<size_t> vuSrcIndices;
			vuSrcIndices.resize(this->vuDimLengths[0]);

			vector<size_t> vuSub;
			// ADD-BY-LEETEN 10/06/2012-END

			/*	// ADD-BY-LEETEN 10/01/2012-BEGIN
			for l = 1 ... log N
				for each dimension, update the first 2^l columns (skip the cell whose subscripts is within [1,2^l]^d)
			*/
			/*
			for(size_t uNrOfRows = 1; uNrOfRows < this->vuDimLengths[0]/2; uNrOfRows <<= 1)
			{
				for(size_t d = 0; d < this->UGetNrOfDims(); d++)
				{
					for(size_t r = 0; r < uNrOfRows; r++)
					{


					}
				}
			}
			*/	// ADD-BY-LEETEN 10/01/2012-BEGIN
			// ADD-BY-LEETEN 10/01/2012-BEGIN
			// DEL-BY-LEETEN 10/08/2012:	vvdBinIsotropicCoefs.resize(this->UGetNrOfBins());
			for(size_t b = 0; b < this->UGetNrOfBins(); b++)
			{
				// DEL-BY-LEETEN 10/08/2012:	vvdBinIsotropicCoefs[b] = this->vvdBinCoefs[b];
				
				// update the hyperslice
				for(size_t uNrOfSlices = 1; uNrOfSlices < this->vuDimLengths[0]/2; uNrOfSlices <<= 1)
				{
					for(size_t d = 0; d < this->UGetNrOfDims(); d++)
					{
						size_t uNrOfElementsInSlice = 1;
						for(size_t d2 = 0; d2 < this->UGetNrOfDims(); d2++)
							uNrOfElementsInSlice *= ( d == d2 )?1:this->vuDimLengths[d2];

						for(size_t e = 0; e < uNrOfElementsInSlice; e++)
						{
							bool bIsValid = false;

							// convert the element to its subscripts
							for(size_t d2 = 0, uIndex = e; d2 < this->UGetNrOfDims(); d2++)
							{
								if( d == d2 )
								{
									vuBase[d2] = 0;	
								}
								else
								{
									size_t uSub = uIndex % this->vuDimLengths[d2];
									vuBase[d2] = uSub;
									uIndex /= this->vuDimLengths[d2];
									if( uSub >= uNrOfSlices * 2 )
										bIsValid = true;
								}
							} // for d2

							if( !bIsValid )
								continue;

							for(size_t s = 0; s < 2 * uNrOfSlices ; s++)
							{
								vuSrc = vuBase; 
								vuSrc[d] = s;
								size_t uSrc = this->UConvetSubToIndex(vuSrc); // MOD-BY-LEETEN 10/05/2012: size_t uSrc = UConvetSubToIndex(vuSrc);
								vdSrc[s] = this->DGetBinCoef(b, uSrc);	// MOD-BY-LEETEN 10/08/2012-FROM:	vdSrc[s] = vvdBinIsotropicCoefs[b][uSrc];
								vuSrcIndices[s] = uSrc;
							}

							for(size_t sp = 0, s = 0; s < uNrOfSlices ; s++)
								for(size_t s2 = 0; s2 < 2; s2++, sp++)
									#if	0	// MOD-BY-LEETEN 10/08/2012-FROM:
										vvdBinIsotropicCoefs[b][vuSrcIndices[sp]] = 
										(	vdSrc[s] + 
											((!s2)?(+1.0):(-1.0)) * vdSrc[s + uNrOfSlices])/2.0;
									#else		// MOD-BY-LEETEN 10/08/2012-TO:
									this->_SetBinCoef(
										b, 
										vuSrcIndices[sp],
										(vdSrc[s] +((!s2)?(+1.0):(-1.0)) * vdSrc[s + uNrOfSlices])/2.0);
									#endif		// MOD-BY-LEETEN 10/08/2012-END
						} // for e
					} // for d
				} // for uNrOfSlices
				
				// ADD-BY-LEETEN 10/05/2012-BEGIN
				// now apply wavelet
#if 0 // MOD-BY-LEETEN 10/09/2012-FROM:
				for(size_t e = 0; e < vvdBinCoefs[b].size(); e++)	// MOD-BY-LEETEN 10/08/2012-FROM:	for(size_t e = 0; e < vvdBinIsotropicCoefs[b].size(); e++)
#else // MOD-BY-LEETEN 10/09/2012-TO:
				for(size_t e = 0; e < this->vvdBinCoefs[b].size(); e++)	
#endif // MOD-BY-LEETEN 10/09/2012-END
				{
					this->_ConvetIndexToSub(e, vuSub);
					size_t uMaxSub = 0;
					for(size_t d = 0; d < vuSub.size(); d++)
						uMaxSub = max(uMaxSub, vuSub[d]);
					double dWaveletWeight = 1.0; 
					if( uMaxSub )
					{
						size_t uLevel = (size_t)ceil(log( (double)(uMaxSub + 1) ) / log(2.0) );
						dWaveletWeight = pow( (double)(1 << (uLevel - 1)), (double)this->UGetNrOfDims()/2.0);
					}
#if 0 // MOD-BY-LEETEN 10/09/2012-FROM:
					vvdBinCoefs[b][e] *= dWaveletWeight / this->dWaveletDenomiator;	// MOD-BY-LEETEN 10/08/2012-FROM:	vvdBinIsotropicCoefs[b][e] *= dWaveletWeight / this->dWaveletDenomiator;
#else // MOD-BY-LEETEN 10/09/2012-TO:
					this->vvdBinCoefs[b][e] *= dWaveletWeight / this->dWaveletDenomiator;	
#endif // MOD-BY-LEETEN 10/09/2012-END
				} // for e

				// ADD-BY-LEETEN 10/08/2012-BEGIN
#if 0 // MOD-BY-LEETEN 10/09/2012-FROM:
				for(map<size_t, double>::iterator
					ipairCoef = vmapBinCoefs[b].begin();
					ipairCoef != vmapBinCoefs[b].end();
					ipairCoef++)
#else // MOD-BY-LEETEN 10/09/2012-TO:
				for(map<size_t, double>::iterator
					ipairCoef = this->vmapBinCoefs[b].begin();
					ipairCoef != this->vmapBinCoefs[b].end();
					ipairCoef++)
#endif // MOD-BY-LEETEN 10/09/2012-END
				{
					size_t e = ipairCoef->first;

					this->_ConvetIndexToSub(e, vuSub);
					size_t uMaxSub = 0;
					for(size_t d = 0; d < vuSub.size(); d++)
						uMaxSub = max(uMaxSub, vuSub[d]);
					double dWaveletWeight = 1.0; 
					if( uMaxSub )
					{
						size_t uLevel = (size_t)ceil(log( (double)(uMaxSub + 1) ) / log(2.0) );
						dWaveletWeight = pow( (double)(1 << (uLevel - 1)), (double)this->UGetNrOfDims()/2.0);
					}
					ipairCoef->second *= dWaveletWeight / this->dWaveletDenomiator;
				}
				// ADD-BY-LEETEN 10/08/2012-END
				// ADD-BY-LEETEN 10/05/2012-END
			} // for b
			// ADD-BY-LEETEN 10/01/2012-END
		}
	};
}
