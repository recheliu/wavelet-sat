#pragma once

#include <valarray>	// ADD-BY-LEETEN 12/30/2012

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
		typename vector<T>::iterator ivSrcHalf = vSrc.begin();
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

	// ADD-BY-LEETEN 12/29/2012-BEGIN
	template<typename T>
	void
	_IDWT1D
	(
		valarray<T>& vSrc,
		valarray<T>& vDst,
		size_t uLength,
		size_t uLevel,
		void* _Reserved = NULL
	)
	{
		size_t uHalfLength = uLength / 2;
		for(size_t c = 0; c < uHalfLength; c++)
		{
			T Src1 = vSrc[c];
			T Src2 = vSrc[c + uHalfLength];
			vDst[c * 2]		= (Src1	+ Src2) * M_SQRT1_2;
			vDst[c * 2 + 1]	= (Src1	- Src2) * M_SQRT1_2;
		}

		if( uLevel > 0 )
		{
			for(size_t c = 0; c < uLength; c++)
				vSrc[c] = vDst[c];
			_IDWT1D<T>(vSrc, vDst, uLength * 2, uLevel - 1);
		}
	}	
	// ADD-BY-LEETEN 12/29/2012-END

	// ADD-BY-LEETEN 01/21/2013-BEGIN
	template<typename T>
	void
	_IDWT1D
	(
		valarray<T>& vSrc,
		valarray<T>& vTemp,
		valarray<T>& vDst,
		size_t uWaveletLength,
		size_t uDataLength,
		size_t uLength,
		size_t uLevel,
		void* _Reserved = NULL
	)
	{
		size_t uHalfLength = uLength / 2;
		size_t uBound = min(uHalfLength, (size_t)ceil((double)uDataLength/(double)uWaveletLength));
		for(size_t c = 0; c < uBound; c++)
		{
			T Src1 = vTemp[c];
			T Src2 = vSrc[c + uHalfLength];
			vDst[c * 2]		= (Src1	+ Src2) * M_SQRT1_2;
			vDst[c * 2 + 1]	= (Src1	- Src2) * M_SQRT1_2;
		}

		if( uLevel > 0 )
		{
			_IDWT1D<T>(
				vSrc, 
				vDst, 
				vTemp,
				uWaveletLength / 2,
				uDataLength,
				uLength * 2, 
				uLevel - 1);
		}
	}	

	// ADD-BY-LEETEN 12/29/2012-BEGIN
	template<typename T>
	void
	_IDWT1D
	(
		const T* pSrc,
		size_t uStep,
		T* pTemp,
		T* pDst,
		size_t uWaveletLength,
		size_t uDataLength,
		size_t uLength,
		size_t uLevel,
		void* _Reserved = NULL
	)
	{
		size_t uHalfLength = uLength / 2;
		// uBound: The bound covered by the data domain.
		size_t uBound = min(uHalfLength, (size_t)ceil((double)uDataLength/(double)uWaveletLength));
		for(size_t c = 0, s = uHalfLength*uStep; c < uBound; c++, s+= uStep)
		{
			T Src1 = pTemp[c];
			T Src2 = pSrc[s];
			pDst[c * 2]		= (Src1	+ Src2) * M_SQRT1_2;
			pDst[c * 2 + 1]	= (Src1	- Src2) * M_SQRT1_2;
		}

		if( uLevel > 0 )
		{
			_IDWT1D<T>(
				pSrc, 
				uStep,
				pDst, // now swap the dst and the temp arrays
				pTemp,
				uWaveletLength / 2,
				uDataLength,
				uLength * 2, 
				uLevel - 1);
		}
	}	
	// ADD-BY-LEETEN 01/21/2013-END
}
