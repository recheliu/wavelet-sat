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
		  // MOD-BY-LEETEN 10/02/2012-FROM:		CBase::_GetAllSums
		  // MOD-BY-LEETEN 10/02/2012-TO:
		  CBase<T>::_GetAllSums
		  // MOD-BY-LEETEN 10/02/2012-END
			(
				vuPos, 
				vdSums
			);
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
							vector<size_t> vuBase;
							for(size_t d2 = 0, uIndex = e; d2 < this->UGetNrOfDims(); d2++)
							{
								if( d == d2 )
								{
									vuBase.push_back(0);
								}
								else
								{
									size_t uSub = uIndex % this->vuDimLengths[d2];
									vuBase.push_back(uSub);
									uIndex /= this->vuDimLengths[d2];
									if( uSub >= uNrOfSlices * 2 )
										bIsValid = true;
								}
							} // for d2

							if( !bIsValid )
								continue;

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
						} // for e
					} // for d
				} // for uNrOfSlices
			} // for b
			// ADD-BY-LEETEN 10/01/2012-END
		}
	};
}
