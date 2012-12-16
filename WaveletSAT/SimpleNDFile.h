#pragma once

#include "SATEncoder.h"
#include "WaveletSATEncoder.h"
#include "SATSepDWTFile.h"
#if WITH_NETCDF // ADD-BY-LEETEN 10/29/2012	
#include <netcdf.h>
#endif

template<class T>
class CSimpleNDFile:
	public WaveletSAT::CSATSepDWTFile<T>
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
};

