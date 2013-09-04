#pragma once

#include "SATSepDWTDecoder.h"
#include "SATFileDecoder.h"	// ADD-BY-LEETEN 01/02/2013
#if WITH_NETCDF // ADD-BY-LEETEN 10/29/2012	
#include <netcdf.h>
#endif

// ADD-BY-LEETEN 12/30/2012-BEGIN
//! Scan the histogra bin to update the entropy
template<class DT>	//!< Type of the data items
DT
ScanHistogramForEntropy
(
	const DT& v
)
{
	return ( v <= (DT)0 )?(DT)0:(DT)(v * log((double)v));
}

//! Clamp the negative value
template<class DT>	//!< Type of the data items
DT
_ClampNegative
(
	const DT& v
)
{
	return max(v, (DT)0);
}
// ADD-BY-LEETEN 12/30/2012-END
template<
	typename DT,							//!< Type of the data
	typename ST = WaveletSAT::typeSum,		//!< Type of the sum
	typename BT = WaveletSAT::typeBin,		//!< Type of the bin
	typename WT = WaveletSAT::typeWavelet	//!< Type of the wavelet coefficientsd
>
class CSimpleNDFile:
	#if	!WITH_SAT_FILE	// ADD-BY-LEETEN 01/05/2013
	virtual public WaveletSAT::CSATSepDWTDecoder<ST, BT, WT>
	// ADD-BY-LEETEN 01/05/2013-BEGIN
	#else	// #if	!WITH_SAT_FILE	
	virtual public WaveletSAT::CSATFileDecoder<ST, BT>
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
	  const size_t& uNrOfBins = this->uNrOfBins; // ADD-BY-LEETEN 01/03/2013

		DT clampedValue = min(max(value, valueMin), valueMax);
		size_t uBin = (size_t)floor((double)((double)uNrOfBins * (double)(clampedValue - valueMin))/(double)(valueMax - valueMin));
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

	// ADD-BY-LEETEN 12/29/2012-BEGIN
	void
	_DecodeSAT(
		void *_Reserved = NULL
	)
	{
		fill(this->vusCachedNextOffsets.begin(), this->vusCachedNextOffsets.end(), (BT)0); 
		// fill(vusCachedBins.begin(),	vusCachedBins.end(), (unsigned short)0);
		// fill(vCachedValues.begin(),	vCachedValues.end(), (WT)0);
		vector<ST> vSAT;
		for(size_t b = 0; b < this->UGetNrOfBins(); b++)
		{
			// LOG_VAR(b);
			_DecodeBin((WaveletSAT::typeBin)b, vSAT);
		}
	}
	// ADD-BY-LEETEN 12/29/2012-END

	// ADD-BY-LEETEN 12/30/2012-BEGIN
	void
	_ComputeEntropy(
		vector<int> viLeft,
		vector<int> viRight,
		// ADD-BY-LEETEN 12/30/2012-BEGIN
		vector<ST>& vEntropyField,
		// ADD-BY-LEETEN 12/30/2012-END
		void *_Reserved = NULL
	)
	{
		// ADD-BY-LEETEN 01/02/2013-BEGIN
		LIBCLOCK_INIT(this->bIsPrintingDecodeBinTiming, __FUNCTION__);
		LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
		// ADD-BY-LEETEN 01/02/2013-END
	  // ADD-BY-LEETEN 12/30/2012-BEGIN
		const size_t uDataSize = this->uDataSize;
		vector<size_t> vuDecodedLengths;

		this->_GetDecodedSize(vuDecodedLengths);

		// ADD-BY-LEETEN 12/30/2012-END
		////////////////////////////////////////////////////////
		// decide the offsets
		size_t uNrOfCorners = (size_t)1<<this->UGetNrOfDims();

		vector< long long > vllOffsets;	vllOffsets.resize(uNrOfCorners);
		vector< int > viSigns;			viSigns.resize(uNrOfCorners);
		for(size_t i = 0; i < uNrOfCorners; i++)
		{
			int iSign = 1; 
			long long llOffset = 0;
			for(size_t 
				d = this->UGetNrOfDims(), j = i; 
				d > 0 ; 
				d--, j /= 2)
			{
				bool bIsLower = (0 == j % 2)?0:1;
				iSign *= (bIsLower)?(-1):(+1);
				llOffset = (llOffset * (long long)vuDecodedLengths[d - 1]) + (long long)((bIsLower)?viLeft[d - 1]:viRight[d - 1]);
			}
			viSigns[i] = iSign;
			vllOffsets[i] = llOffset;
		}
		LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);	// ADD-BY-LEETEN 01/02/2013

		/////////////// compute the SAT
		LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);	// ADD-BY-LEETEN 01/02/2013
		vector<ST> vSAT;	// ADD-BY-LEETEN 01/09/2013
		vector<ST> vLocalHist;	
		vector<ST> vTempEntropyField;	
		vector<ST> vSum;			
		for(size_t b = 0; b < this->UGetNrOfBins(); b++)
		{
			// ADD-BY-LEETEN 01/02/2013-BEGIN
			LIBCLOCK_INIT(this->bIsPrintingDecodeBinTiming, __FUNCTION__);
			LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
			// ADD-BY-LEETEN 01/02/2013-END

			_DecodeBin((BT)b, vSAT);

			// ADD-BY-LEETEN 01/02/2013-BEGIN
			LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);

			LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
			if( vSAT.size() != vLocalHist.size() )
				vLocalHist.resize(vSAT.size());
			if( vSAT.size() != vTempEntropyField.size() )
				vTempEntropyField.resize(vSAT.size());
			if( vSAT.size() != vSum.size() )
				vSum.resize(vSAT.size());
			LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);
			// ADD-BY-LEETEN 01/02/2013-END

			// compute the local sum
			vLocalHist.assign(vLocalHist.size(), (ST)0);

			LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);	// ADD-BY-LEETEN 01/02/2013
			for(size_t i = 0; i < uNrOfCorners; i++)
			{
				size_t uBegin, uEnd;
				if( !vllOffsets[i] )
					continue;
				if( vllOffsets[i] > 0 )
				{
					uBegin = 0;
					uEnd = vSAT.size() - vllOffsets[i]; 
				}
				if( vllOffsets[i] < 0 )
				{
					uBegin = (size_t)-vllOffsets[i];
					uEnd = vSAT.size(); 
				}
				for(size_t d = uBegin; d < uEnd; d++)
					vLocalHist[d] += vSAT[d + vllOffsets[i]] * (ST)viSigns[i];
			}
			LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);

			LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
			for(size_t d = 0; d < vLocalHist.size(); d++)
			{
				if( vLocalHist[d] <= (DT)0 )
					continue;
				vTempEntropyField[d] += ScanHistogramForEntropy(vLocalHist[d]);
				vSum[d] += vLocalHist[d];
			}

			// ADD-BY-LEETEN 01/02/2013-BEGIN
			LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);
			LIBCLOCK_PRINT(this->bIsPrintingDecodeBinTiming);
			// ADD-BY-LEETEN 01/02/2013-END
		}
		for(size_t i = 0; i < vTempEntropyField.size(); i++)
		{
			vTempEntropyField[i] = -vTempEntropyField[i] / vSum[i] + log(vSum[i]);
			vTempEntropyField[i] /= (ST)M_LN2;
		}
		LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);	// ADD-BY-LEETEN 01/02/2013

		LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
		_ClampToDataSize(vTempEntropyField, vEntropyField);
		_ClampBorder(vEntropyField, viLeft, viRight);
		LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);
		LIBCLOCK_PRINT(this->bIsPrintingDecodeBinTiming);
	}
	// ADD-BY-LEETEN 12/30/2012-END
	// ADD-BY-LEETEN 2013/09/03-BEGIN
	#if	!WITH_SAT_FILE
	// Given the level and local coordinate within the level, 
	// return the wavelet coefficients
	void 
	_GetCoefSparse
	(
		const vector<size_t>& vuLevel,
		const vector<size_t>& vuLocal,
		vector< pair<BT, WT> >& vpairCoefBinValues,
		void* _Reserved = NULL
	)
	{
		vpcCoefPools.at(WaveletSAT::UConvertSubToIndex(vuLevel, vuDimLevels))->_GetCoefSparse(
			vuLocal, 
			vpairCoefBinValues);
	}
	#endif	// #if	!WITH_SAT_FILE
	// ADD-BY-LEETEN 2013/09/03-END
};

