#pragma once

#pragma warning (disable:	4250	4996)

#include <vector>
using namespace std;
#include <math.h>

#include "utils.h" 

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	typedef unsigned short	typeBin;
	#if	WITH_DOUBLE_COEF
	typedef double			typeWavelet;
	#else	// #if	WITH_DOUBLE_COEF
	typedef float	typeWavelet;
	#endif	// #if	WITH_DOUBLE_COEF
	typedef double			typeSum;

	//! The base class of all classes
	class CBase
	{
protected:	
public:
		enum EParameter
		{
			PARAMETER_BEGIN = 0x0000,
			PARAMETER_END
		};

		virtual
		void
		_SetInteger(
			int eName,
			long lValue,
			void* _Reserved = NULL
		)
		{
		}

		virtual
		void
		_GetInteger(
			int eName,
			long *lValue,
			void* _Reserved = NULL
		)
		{
		}
	};
}

/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
