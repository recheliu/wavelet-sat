#pragma once

#include "SATEncoder.h"
#include "WaveletSATEncoder.h"
#include "SATSepDWTEncoder.h"
#include "SATFileEncoder.h"	// ADD-BY-LEETEN 10/28/2012

template<class T>
class CSimpleND:
// ADD-BY-LEETEN 10/28/2012-BEGIN
#if	WITH_NETCDF
	public WaveletSAT::CSATFileEncoder<T, double>
#else // #if	WITH_NETCDF
// ADD-BY-LEETEN 10/28/2012-END
public WaveletSAT::CWaveletSATEncoder<T>
#endif	// #if	WITH_NETCDF	// ADD-BY-LEETEN 10/28/2012
	// MOD-BY-LEETEN 10/28/2012-END
{
	size_t uNrOfBins;
	T valueMin, valueMax;
public:
	virtual 
	void 
	_MapValueToBins
	(
		const vector<size_t>& vuPos,
		const T& value, 
		vector< pair<size_t, double> >& vpBins,
		void *_Reserved = NULL
	)
	{
		T clampedValue = min(max(value, valueMin), valueMax);
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
		const T& valueMin, 
		const T& valueMax
	)
	{
		this->uNrOfBins = uNrOfBins;
		this->valueMin = valueMin;
		this->valueMax = valueMax;
	}

	void
	_AddValue
	(
		const vector<size_t>& vuPos,
		const T& value,
		void *_Reserved = NULL
	)
	{
		_Update(vuPos, value);
	}
};

