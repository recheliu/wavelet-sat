#pragma once

#if	0	// DEL-BY-LEETEN 12/25/2012-BEGIN
#include "SATEncoder.h"
#include "WaveletSATEncoder.h"
#endif	// DEL-BY-LEETEN 12/25/2012-END
// MOD-BY-LEETEN 12/25/2012-FROM:	#include "SATSepDWTFile.h"
#include "SATSepDWTOutOfCore.h"
// MOD-BY-LEETEN 12/25/2012-END
#if WITH_NETCDF // ADD-BY-LEETEN 10/29/2012	
#include <netcdf.h>
#endif

// MOD-BY-LEETEN 12/29/2012-FROM:	template<class DT>
template<
	class DT,	//!< Type of the data items
	class WT	//!< Type of the wavelet coefficients.
>
// MOD-BY-LEETEN 12/29/2012-END
class CSimpleNDFile:
	// MOD-BY-LEETEN 12/25/2012-FROM:	public WaveletSAT::CSATSepDWTFile<DT>
	// MOD-BY-LEETEN 12/29/2012-FROM:	public WaveletSAT::CSATSepDWTOutOfCore<DT>
	public WaveletSAT::CSATSepDWTOutOfCore<DT, WT>
	// MOD-BY-LEETEN 12/29/2012-END
	// MOD-BY-LEETEN 12/25/2012-END
{
	size_t uNrOfBins;
	DT valueMin, valueMax;
public:
	virtual 
	void 
	_MapValueToBins
	(
		const vector<size_t>& vuPos,
		const DT& value, 
		vector< pair<size_t, double> >& vpBins,
		void *_Reserved = NULL
	)
	{
		DT clampedValue = min(max(value, valueMin), valueMax);
		size_t uBin = (size_t)floorf((float)(uNrOfBins * (clampedValue - valueMin))/(float)(valueMax - valueMin));
		uBin = min(uBin, uNrOfBins - 1);
		vpBins.clear();
		vpBins.push_back(pair<size_t, double>(uBin, 1.0));
	}

	////////////////////////////////////////////////////////////
	void
	_SetHistogram
	(
		size_t uNrOfBins,
		const DT& valueMin, 
		const DT& valueMax
	)
	{
		this->uNrOfBins = uNrOfBins;
		this->valueMin = valueMin;
		this->valueMax = valueMax;
	}
};

