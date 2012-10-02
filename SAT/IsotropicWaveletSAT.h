#pragma once

// ADD-BY-LEETEN 09/12/2012-BEGIN
// DEL-BY-LEETEN 09/13:	#define WITH_SPARSE_WAVELET_COEFS	0

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
			enum EParameter eName,
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
			CBase::_SetLong(eName, lValue);
			// ADD-BY-LEETEN 10/01/2012-END
		}
		
		// ADD-BY-LEETEN 10/01/2012-BEGIN
		virtual	
		void
		_ShowStatistics
		(
		)
		{
			CBase::_ShowStatistics();

			size_t uNrOfNonZeroCoefs = 0;
			for(size_t b = 0; b < UGetNrOfBins(); b++)
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
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				LOG_VAR(vuDimLengths[d]);
				uNrOfDataItems *= this->vuDataDimLengths[d];	// MOD-BY-LEETEN 09/30/2012-FROM:	uNrOfDataItems *= this->vuDimLengths[d];
			}
			LOG_VAR(uNrOfDataItems);
			LOG_VAR(UGetNrOfBins());

			double dCR = (double)(uNrOfDataItems * UGetNrOfBins()) / (double)uNrOfNonZeroCoefs;
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
			CBase::_GetAllSums
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
			// ADD-BY-LEETEN 10/01/2012-BEGIN
			// Do not apply wavelet coefficients
			CBase::_SetBoolean(CBase::FINALIZED_WITHOUT_WAVELET, true);
			// ADD-BY-LEETEN 10/01/2012-END
			CBase::_Finalize
			(
			);

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
				vvdBinIsotropicCoefs[b] = vvdBinCoefs[b];
				
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
