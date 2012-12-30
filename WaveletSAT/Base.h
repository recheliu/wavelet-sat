#pragma once

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
	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

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
