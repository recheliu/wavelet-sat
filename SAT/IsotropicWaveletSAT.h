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
		vector< vector<double> > vvdBinIsotropicCoefs;

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
			// MOD-BY-LEETEN 10/02/2012-FROM:			enum EParameter eName,
			// MOD-BY-LEETEN 10/02/2012-TO:
			int eName,
			// MOD-BY-LEETEN 10/02/2012-END
			long lValue,
			void* _Reserved = NULL
		)
		{
			// ADD-BY-LEETEN 10/01/2012-BEGIN
			switch(eName)
			{
			default:
				;
			}
			// MOD-BY-LEETEN 10/02/2012-FROM:			CBase::_SetLong(eName, lValue);
			// MOD-BY-LEETEN 10/02/2012-TO:
			CBase<T>::_SetLong(eName, lValue);
			// MOD-BY-LEETEN 10/02/2012-END
			// ADD-BY-LEETEN 10/01/2012-END
		}
		
		// ADD-BY-LEETEN 10/01/2012-BEGIN
		virtual	
		void
		_ShowStatistics
		(
		)
		{
		  // MOD-BY-LEETEN 10/02/2012-FROM:			CBase::_ShowStatistics();
		  // MOD-BY-LEETEN 10/02/2012-TO:
		  CBase<T>::_ShowStatistics();
		  // MOD-BY-LEETEN 10/02/2012-END

			size_t uNrOfNonZeroCoefs = 0;
			// MOD-BY-LEETEN 10/02/2012-FROM:			for(size_t b = 0; b < UGetNrOfBins(); b++)
			// MOD-BY-LEETEN 10/02/2012-TO:
			for(size_t b = 0; b < this->UGetNrOfBins(); b++)
			// MOD-BY-LEETEN 10/02/2012-END
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
			// MOD-BY-LEETEN 10/02/2012-FROM:			for(size_t d = 0; d < UGetNrOfDims(); d++)
			// MOD-BY-LEETEN 10/02/2012-TO:
			for(size_t d = 0; d < this->UGetNrOfDims(); d++)
			// MOD-BY-LEETEN 10/02/2012-END
			{
#if 0 			// MOD-BY-LEETEN 10/02/2012-FROM:
				LOG_VAR(vuDimLengths[d]);
				uNrOfDataItems *= this->vuDataDimLengths[d];	// MOD-BY-LEETEN 09/30/2012-FROM:	uNrOfDataItems *= this->vuDimLengths[d];
#else			// MOD-BY-LEETEN 10/02/2012-FROM:
				LOG_VAR(this->vuDimLengths[d]);
				uNrOfDataItems *= this->vuDataDimLengths[d];
#endif			// MOD-BY-LEETEN 10/02/2012-END
			}
			LOG_VAR(uNrOfDataItems);
			// MOD-BY-LEETEN 10/02/2012-FROM:			LOG_VAR(UGetNrOfBins());
			// MOD-BY-LEETEN 10/02/2012-TO:
			LOG_VAR(this->UGetNrOfBins());
			// MOD-BY-LEETEN 10/02/2012-END

			// MOD-BY-LEETEN 10/02/2012-FROM:			double dCR = (double)(uNrOfDataItems * UGetNrOfBins()) / (double)uNrOfNonZeroCoefs;
			// MOD-BY-LEETEN 10/02/2012-TO:
			double dCR = (double)(uNrOfDataItems * this->UGetNrOfBins()) / (double)uNrOfNonZeroCoefs;
			// MOD-BY-LEETEN 10/02/2012-END
			LOG_VAR(dCR);

			double dOverhead = (double)uNrOfNonZeroCoefs / (double)uNrOfDataItems;
			LOG_VAR(dOverhead);
		}
		// ADD-BY-LEETEN 10/01/2012-END

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
		#if	0	// MOD-BY-LEETEN 10/05/2012-FROM:
			  // MOD-BY-LEETEN 10/02/2012-FROM:		CBase::_GetAllSums
			  // MOD-BY-LEETEN 10/02/2012-TO:
			  CBase<T>::_GetAllSums
			  // MOD-BY-LEETEN 10/02/2012-END
				(
					vuPos, 
					vdSums
				);
		#else		// MOD-BY-LEETEN 10/05/2012-TO:
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
				dCount += this->vvdBinIsotropicCoefs[b][0];

				// now get the wavelet coefficients
				for(size_t 
					l = 0, uWaveletSize = (size_t)1 << (this->vuDimLevels[0] - 1);
					l < (this->vuDimLevels[0] - 1); 
					l++, uWaveletSize >>= 1)
				{
					double dWaveletWeight = pow( (double)(1 << l), (double)this->UGetNrOfDims()/2.0);

					// get the wavelet position
					#if	0	// DEL-BY-LEETEN 10/06/2012-BEGIN
						vector<size_t> vuWaveletPos;
						vector<size_t> vuPosInWavelet;
					#endif		// DEL-BY-LEETEN 10/06/2012-END
					for(size_t d = 0; d < vuPos.size(); d++)
					{
						#if	0	// MOD-BY-LEETEN 10/06/2012-FROM:
							vuWaveletPos.push_back( vuPos[d] / uWaveletSize );
							vuPosInWavelet.push_back( vuPos[d] % uWaveletSize );
						#else		// MOD-BY-LEETEN 10/06/2012-TO:
						vuWaveletPos[d] = vuPos[d] / uWaveletSize;
						vuPosInWavelet[d] = vuPos[d] % uWaveletSize;
						#endif		// MOD-BY-LEETEN 10/06/2012-END
					}

					double dWaveletMaxAbsCoef = 0.0;	// ADD-BY-LEETEN 10/06/2012

					for(size_t p = 0, w = 0; w < uNrOfWaveletsPerLevel; w++, p += this->UGetNrOfDims())
					{
						// locate the wavelet coefficients
						#if	0	// MOD-BY-LEETEN 10/06/2012-FROM:
							vector< size_t > vuCoefPos;
							for(size_t d = 0; d < this->UGetNrOfDims(); d++)
								vuCoefPos.push_back(vuWaveletToOffset[p + d] * ((size_t) 1 << l) + vuWaveletPos[d]);
						#else		// MOD-BY-LEETEN 10/06/2012-TO:
						for(size_t d = 0; d < this->UGetNrOfDims(); d++)
							vuCoefPos[d] = vuWaveletToOffset[p + d] * ((size_t) 1 << l) + vuWaveletPos[d];
						#endif		// MOD-BY-LEETEN 10/06/2012-END

						size_t uIndex = this->UConvetSubToIndex(vuCoefPos);
						double dWaveletCoef = this->vvdBinIsotropicCoefs[b][uIndex];

						// now find the basis
#if 0 // MOD-BY-LEETEN 10/06/2012-FROM:
						if( fabs(dWaveletCoef) >= dWaveletThreshold )	// MOD-BY-LEETEN 10/06/2012-FROM:	if( dWaveletCoef )
#else // MOD-BY-LEETEN 10/06/2012-TO:
						if( fabs(dWaveletCoef) >= this->dWaveletThreshold )	// MOD-BY-LEETEN 10/06/2012-FROM:	if( dWaveletCoef )
#endif // MOD-BY-LEETEN 10/06/2012-END
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
		#endif		// MOD-BY-LEETEN 10/05/2012-END
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
#if 0 		// MOD-BY-LEETEN 10/02/2012-FROM:
			// ADD-BY-LEETEN 10/01/2012-BEGIN
			// Do not apply wavelet coefficients
			CBase::_SetBoolean(CBase::FINALIZED_WITHOUT_WAVELET, true);
			// ADD-BY-LEETEN 10/01/2012-END
			CBase::_Finalize
			(
			);
#else		// MOD-BY-LEETEN 10/02/2012-TO:
			// Do not apply wavelet coefficients
			this->_SetBoolean(CBase<T>::FINALIZED_WITHOUT_WAVELET, true);
			CBase<T>::_Finalize
			(
			);
#endif		// MOD-BY-LEETEN 10/02/2012-END

			// ADD-BY-LEETEN 10/06/2012-BEGIN
			vector<size_t> vuBase;
			vuBase.resize(this->UGetNrOfDims()); // MOD-BY-LEETEN 10/06/2012-FROM: .resize(UGetNrOfDims());

			vector<size_t> vuSrc;
			vuSrc.resize(this->UGetNrOfDims()); // MOD-BY-LEETEN 10/06/2012-FROM: vuSrc.resize(UGetNrOfDims());

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
			vvdBinIsotropicCoefs.resize(this->UGetNrOfBins());
			for(size_t b = 0; b < this->UGetNrOfBins(); b++)
			{
			  // MOD-BY-LEETEN 10/02/2012-FROM:				vvdBinIsotropicCoefs[b] = vvdBinCoefs[b];
			  // MOD-BY-LEETEN 10/02/2012-TO:
			  vvdBinIsotropicCoefs[b] = this->vvdBinCoefs[b];
			  // MOD-BY-LEETEN 10/02/2012-END
				
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
							// DEL-BY-LEETEN 10/06/2012:	vector<size_t> vuBase;
							for(size_t d2 = 0, uIndex = e; d2 < this->UGetNrOfDims(); d2++)
							{
								if( d == d2 )
								{
									vuBase[d2] = 0;	// MOD-BY-LEETEN 10/06/2012-FROM:	vuBase.push_back(0);
								}
								else
								{
									size_t uSub = uIndex % this->vuDimLengths[d2];
									vuBase[d2] = uSub;	// MOD-BY-LEETEN 10/06/2012-FROM:	vuBase.push_back(uSub);
									uIndex /= this->vuDimLengths[d2];
									if( uSub >= uNrOfSlices * 2 )
										bIsValid = true;
								}
							} // for d2

							if( !bIsValid )
								continue;

							#if	0	// MOD-BY-LEETEN 10/05/2012-FROM:
								for(size_t s = 0; s < uNrOfSlices ; s++)
								{
	#if 0 							// MOD-BY-LEETEN 10/02/2012-FROM:
									vector<size_t> vuSrc1  = vuBase; vuSrc1[d] = s;
									size_t uSrc1 = UConvetSubToIndex(vuSrc1);
									double dSrc1 = vvdBinCoefs[b][uSrc1];

									vector<size_t> vuSrc2  = vuBase; vuSrc2[d] = s + uNrOfSlices;
									size_t uSrc2 = UConvetSubToIndex(vuSrc2);
									double dSrc2 = vvdBinCoefs[b][uSrc2];

									vector<size_t> vuDst1  = vuBase; vuDst1[d] = s;
									size_t uDst1 = UConvetSubToIndex(vuDst1);

									vector<size_t> vuDst2  = vuBase; vuDst2[d] = s + 1;
									size_t uDst2 = UConvetSubToIndex(vuDst2);
	#else							// MOD-BY-LEETEN 10/02/2012-TO:
									vector<size_t> vuSrc1  = vuBase; vuSrc1[d] = s;
									size_t uSrc1 = this->UConvetSubToIndex(vuSrc1);
									double dSrc1 = this->vvdBinCoefs[b][uSrc1];

									vector<size_t> vuSrc2  = vuBase; vuSrc2[d] = s + uNrOfSlices;
									size_t uSrc2 = this->UConvetSubToIndex(vuSrc2);
									double dSrc2 = this->vvdBinCoefs[b][uSrc2];

									vector<size_t> vuDst1  = vuBase; vuDst1[d] = s;
									size_t uDst1 = this->UConvetSubToIndex(vuDst1);

									vector<size_t> vuDst2  = vuBase; vuDst2[d] = s + 1;
									size_t uDst2 = this->UConvetSubToIndex(vuDst2);
	#endif							// MOD-BY-LEETEN 10/02/2012-END

									// Dst1 = Src1 + Src2
									vvdBinIsotropicCoefs[b][uDst1] = dSrc1 + dSrc2;

									// Dst2 = Src1 - Src2
									vvdBinIsotropicCoefs[b][uDst2] = dSrc1 - dSrc2;
								} // for s
							#else	// MOD-BY-LEETEN 10/05/2012-TO:
							#if	0	// MOD-BY-LEETEN 10/06/2012-FROM:
							vector< double > vdSrc;
							for(size_t s = 0; s < 2 * uNrOfSlices ; s++)
							{
								vector<size_t> vuSrc  = vuBase; vuSrc[d] = s;
								size_t uSrc = this->UConvetSubToIndex(vuSrc); // MOD-BY-LEETEN 10/05/2012: size_t uSrc = UConvetSubToIndex(vuSrc);
								vdSrc.push_back( vvdBinIsotropicCoefs[b][uSrc] );
							}

							vector< double > vdDst;
							for(size_t s = 0; s < uNrOfSlices ; s++)
							{
								vdDst.push_back( (vdSrc[s] + vdSrc[s + uNrOfSlices])/2.0 );
								vdDst.push_back( (vdSrc[s] - vdSrc[s + uNrOfSlices])/2.0 );
							}

							for(size_t s = 0; s < 2 * uNrOfSlices ; s++)
							{
								vector<size_t> vuDst  = vuBase; vuDst[d] = s;
								size_t uDst = this->UConvetSubToIndex(vuDst); // MOD-BY-LEETEN 10/05/2012: size_t uDst = UConvetSubToIndex(vuDst);
								vvdBinIsotropicCoefs[b][uDst] = vdDst[s];
							}
							#else		// MOD-BY-LEETEN 10/06/2012-TO:
							for(size_t s = 0; s < 2 * uNrOfSlices ; s++)
							{
								vuSrc = vuBase; 
								vuSrc[d] = s;
								size_t uSrc = this->UConvetSubToIndex(vuSrc); // MOD-BY-LEETEN 10/05/2012: size_t uSrc = UConvetSubToIndex(vuSrc);
								vdSrc[s] = vvdBinIsotropicCoefs[b][uSrc];
								vuSrcIndices[s] = uSrc;
							}

							for(size_t sp = 0, s = 0; s < uNrOfSlices ; s++)
								for(size_t s2 = 0; s2 < 2; s2++, sp++)
									vvdBinIsotropicCoefs[b][vuSrcIndices[sp]] = 
									(	vdSrc[s] + 
										((!s2)?(+1.0):(-1.0)) * vdSrc[s + uNrOfSlices])/2.0;
							#endif		// MOD-BY-LEETEN 10/06/2012-END
							#endif	// MOD-BY-LEETEN 10/05/2012-END
						} // for e
					} // for d
				} // for uNrOfSlices
				
				// ADD-BY-LEETEN 10/05/2012-BEGIN
				// now apply wavelet
				for(size_t e = 0; e < vvdBinIsotropicCoefs[b].size(); e++)
				{
					// DEL-BY-LEETEN 10/06/2012:	vector<size_t> vuSub;
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
					vvdBinIsotropicCoefs[b][e] *= dWaveletWeight / this->dWaveletDenomiator;
				} // for e
				// ADD-BY-LEETEN 10/05/2012-END
			} // for b
			// ADD-BY-LEETEN 10/01/2012-END
		}
	};
}
