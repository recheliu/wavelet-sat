#pragma once

// ADD-BY-LEETEN 12/31/2012-BEGIN
#pragma warning (disable:	4250	4996)
// ADD-BY-LEETEN 12/31/2012-END

#include <vector>
using namespace std;
#include <math.h>

#include "utils.h" // ADD-BY-LEETEN 10/30/2012

#if	0	// DEL-BY-LEETEN 01/18/2012-BEGIN
// ADD-BY-LEETEN 01/05/2013-BEGIN
// If this is non-0, the non-wavelet SAT file will be used
#define WITH_SAT_FILE			0
// ADD-BY-LEETEN 01/05/2013-END
#endif	// DEL-BY-LEETEN 01/18/2012-END

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	// ADD-BY-LEETEN 01/03/2013-BEGIN
	typedef unsigned short	typeBin;
	typedef double			typeWavelet;
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
