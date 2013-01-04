#pragma once

#include "SATEncoder.h"
#include "WaveletSATEncoder.h"
#include "SATSepDWTEncoder.h"
#include "SATFileEncoder.h"		// ADD-BY-LEETEN 01/02/2013
#if WITH_NETCDF // ADD-BY-LEETEN 10/29/2012	
#include <netcdf.h>
#endif

#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
template<class T>
class CSimpleND:
// ADD-BY-LEETEN 10/28/2012-BEGIN
	virtual public WaveletSAT::CWaveletSATEncoder<T>
#else	// MOD-BY-LEETEN 01/03/2013-TO:
template<
	typename DT,							//!< Type of the data
	typename ST = WaveletSAT::typeSum,		//!< Type of the sum
	typename BT = WaveletSAT::typeBin,		//!< Type of the bin
	typename WT = WaveletSAT::typeWavelet	//!< Type of the wavelet coefficientsd
>
class CSimpleND:
	virtual public WaveletSAT::CWaveletSATEncoder<DT, ST, BT, WT>
#endif	// MOD-BY-LEETEN 01/03/2013-END
{
	// MOD-BY-LEETEN 01/03/2013-FROM:	size_t uNrOfBins;
	BT uNrOfBins;
	// MOD-BY-LEETEN 01/03/2013-END
	DT valueMin, valueMax;
public:
	virtual 
	void 
	_MapValueToBins
	(
		const vector<size_t>& vuPos,
		const DT& value, 
		// MOD-BY-LEETEN 01/03/2013-FROM:		vector< pair<size_t, double> >& vpBins,
		vector< pair<BT, ST> >& vpBins,
		// MOD-BY-LEETEN 01/03/2013-END
		void *_Reserved = NULL
	)
	{
		DT clampedValue = min(max(value, valueMin), valueMax);
		#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
		size_t uBin = (size_t)floorf((float)(uNrOfBins * (clampedValue - valueMin))/(float)(valueMax - valueMin));
		uBin = min(uBin, uNrOfBins - 1);
		vpBins.clear();
		vpBins.push_back(pair<size_t, double>(uBin, 1.0));
		#else	// MOD-BY-LEETEN 01/03/2013-TO:
		BT uBin = (BT)floor((double)(uNrOfBins * (clampedValue - valueMin))/(double)(valueMax - valueMin));
		uBin = min(uBin, uNrOfBins - 1);
		vpBins.clear();
		vpBins.push_back(pair<BT, WT>(uBin, (WT)1));
		#endif	// MOD-BY-LEETEN 01/03/2013-END
	}
	
	////////////////////////////////////////////////////////////
	void
	_SetHistogram
	(
		// MOD-BY-LEETEN 01/03/2013-FROM:		size_t uNrOfBins,
		const BT& uNrOfBins,
		// MOD-BY-LEETEN 01/03/2013-END
		const DT& valueMin, 
		const DT& valueMax
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
		const DT& value,
		void *_Reserved = NULL
	)
	{
		_Update(vuPos, value);
	}
};

