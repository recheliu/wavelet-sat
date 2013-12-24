#pragma once

// ADD-BY-LEETEN 12/31/2012-BEGIN
#pragma warning (disable:	4250	4996)
// ADD-BY-LEETEN 12/31/2012-END

#include <vector>
using namespace std;
#include <math.h>

#include "utils.h" // ADD-BY-LEETEN 10/30/2012

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	// ADD-BY-LEETEN 01/03/2013-BEGIN
	typedef unsigned short	typeBin;
	#if	WITH_DOUBLE_COEF
	typedef double			typeWavelet;
	#else	// #if	WITH_DOUBLE_COEF
	typedef float	typeWavelet;
	#endif	// #if	WITH_DOUBLE_COEF
	typedef double			typeSum;
	// ADD-BY-LEETEN 01/03/2013-END
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
