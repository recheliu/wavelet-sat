#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
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
}
