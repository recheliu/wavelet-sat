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
public:
		virtual	// ADD-BY-LEETEN 09/29/2012
		void
		_SetLong(
			enum EParameter eName,
			long lValue,
			void* _Reserved = NULL
		)
		{
		}
		
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
			CBase::_Finalize
			(
			);

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
		}
	};
}
