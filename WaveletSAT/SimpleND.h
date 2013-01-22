#pragma once

#include "WaveletSATEncoder.h"
// ADD-BY-LEETEN 01/10/2012-BEGIN
#if	WITH_CUDA
#include "WaveletSATGPUEncoder.h"	
#endif	// #if	WITH_CUDA
// ADD-BY-LEETEN 01/10/2012-END
#include "SATFileEncoder.h"		// ADD-BY-LEETEN 01/02/2013
#if WITH_NETCDF // ADD-BY-LEETEN 10/29/2012	
#include <netcdf.h>
#endif

template<
	typename DT,							//!< Type of the data
	typename ST = WaveletSAT::typeSum,		//!< Type of the sum
	typename BT = WaveletSAT::typeBin,		//!< Type of the bin
	typename WT = WaveletSAT::typeWavelet	//!< Type of the wavelet coefficientsd
>
class CSimpleND:
	#if	!WITH_SAT_FILE	// ADD-BY-LEETEN 01/05/2013
	#if	!WITH_CUDA		// ADD-BY-LEETEN 01/10/2012
	virtual public WaveletSAT::CWaveletSATEncoder<DT, ST, BT, WT>
	// ADD-BY-LEETEN 01/10/2012-BEGIN
	#else	// #if	!WITH_CUDA	
	virtual public WaveletSAT::CWaveletSATGPUEncoder<DT, ST, BT, WT>
	#endif	// #if	!WITH_CUDA	
	// ADD-BY-LEETEN 01/10/2012-END
	// ADD-BY-LEETEN 01/05/2013-BEGIN
	#else	// #if	!WITH_SAT_FILE	
	virtual public WaveletSAT::CSATFileEncoder<DT, ST, BT>
	#endif	// #if	!WITH_SAT_FILE	
	// ADD-BY-LEETEN 01/05/2013-END
{
	DT valueMin, valueMax;
public:
	virtual 
	void 
	_MapValueToBins
	(
		const vector<size_t>& vuPos,
		const DT& value, 
		vector< pair<BT, ST> >& vpBins,
		void *_Reserved = NULL
	)
	{
	  const size_t& uNrOfBins = this->uNrOfBins; // ADD-BY-LEETEN 01/04/2013

		DT clampedValue = min(max(value, valueMin), valueMax);
		size_t uBin = (size_t)floorf((float)(uNrOfBins * (clampedValue - valueMin))/(float)(valueMax - valueMin));
		uBin = min(uBin, uNrOfBins - 1);
		vpBins.clear();
		vpBins.push_back(pair<BT, ST>((BT)uBin, (ST)1));
	}
	
	////////////////////////////////////////////////////////////
	void
	_SetHistogram
	(
		const DT& valueMin, 
		const DT& valueMax
	)
	{
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

