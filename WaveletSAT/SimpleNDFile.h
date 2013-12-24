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

	// ADD-BY-LEETEN 2013/12/07-BEGIN
	void
	_ComputeRegionSum
	(
		const vector<int>& viLeft,
		const vector<int>& viRight,
		const vector<ST>& vSAT,
		const vector<size_t>& vuLengths,	// ADD-BY-LEETEN 2013/12/14
		vector<ST>& vRegionSum,
		void *_Reserved = NULL
	)
	{
		// ADD-BY-LEETEN 2013/12/14-BEGIN
		if(vRegionSum.size() != vSAT.size() )
			vRegionSum.resize(vSAT.size());
		vRegionSum.assign(vRegionSum.size(), (ST)0);
		// ADD-BY-LEETEN 2013/12/14-END

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
				// MOD-BY-LEETEN 2013/12/14:	llOffset = (llOffset * (long long)vuCoefLengths[d - 1]) + (long long)((bIsLower)?viLeft[d - 1]:viRight[d - 1]);
				llOffset = (llOffset * (long long)vuLengths[d - 1]) + (long long)((bIsLower)?viLeft[d - 1]:viRight[d - 1]);
				// MOD-BY-LEETEN 2013/12/14-END
			}
			viSigns[i] = iSign;
			vllOffsets[i] = llOffset;
		}

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
				vRegionSum[d] += vSAT[d + vllOffsets[i]] * (ST)viSigns[i];
		}
	}
	// ADD-BY-LEETEN 2013/12/07-EN


	// ADD-BY-LEETEN 2013/12/01-BEGIN
	void
	_ComputeEntropy(
		#if	0	// MOD-BY-LEETEN 2013/12/07-FROM:
		vector<int> viLeft,
		vector<int> viRight,
		#else	// MOD-BY-LEETEN 2013/12/07-TO:
		const vector<int>& viLeft,
		const vector<int>& viRight,
		#endif	// MOD-BY-LEETEN 2013/12/07-END
		size_t uNrOfEntropyBins,
		vector<ST>& vEntropyField,
		void *_Reserved = NULL
	)
	{
		LIBCLOCK_INIT(this->bIsPrintingDecodeBinTiming, __FUNCTION__);
		LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
		const size_t uDataSize = this->uDataSize;

		#if	0	// DEL-BY-LEETEN 2013/12/07-BEGIN
		vector<size_t> vuDecodedLengths;

		this->_GetDecodedSize(vuDecodedLengths);

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
		#endif	// DEL-BY-LEETEN 2013/12/07-END
		LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);	

		/////////////// compute the SAT
		LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);	
		vector<ST> vSAT;	
		vector<ST> vLocalHist;	
		vector<ST> vTempEntropyField;	
		vector<ST> vSum;			

		size_t uNrOfAggregatedBins = ( 0 == uNrOfEntropyBins ) ? this->UGetNrOfBins() : uNrOfEntropyBins;
		size_t uBinInterval = this->UGetNrOfBins() / uNrOfAggregatedBins;
		for(size_t b = 0; b < uNrOfAggregatedBins; b++)
		{
			LIBCLOCK_INIT(this->bIsPrintingDecodeBinTiming, __FUNCTION__);
			LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);

			_DecodeAggregatedBin((BT)b * uBinInterval, (BT)(b + 1) * uBinInterval - 1, vSAT);

			LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);

			LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
			if( vSAT.size() != vLocalHist.size() )
				vLocalHist.resize(vSAT.size());
			if( vSAT.size() != vTempEntropyField.size() )
				vTempEntropyField.resize(vSAT.size());
			if( vSAT.size() != vSum.size() )
				vSum.resize(vSAT.size());
			LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);

			// compute the local sum
			vLocalHist.assign(vLocalHist.size(), (ST)0);

			LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
			#if	0	// MOD-BY-LEETEN 2013/12/07-FROM:
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
			#else	// MOD-BY-LEETEN 2013/12/07-TO:
			_ComputeRegionSum
			(
				viLeft,
				viRight,
				vSAT,
				vuCoefLengths,			// ADD-BY-LEETEN 2013/12/14
				vLocalHist
			);
			#endif	// MOD-BY-LEETEN 2013/12/07-END
			LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);

			LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
			for(size_t d = 0; d < vLocalHist.size(); d++)
			{
				if( vLocalHist[d] <= (DT)0 )
					continue;
				vTempEntropyField[d] += ScanHistogramForEntropy(vLocalHist[d]);
				vSum[d] += vLocalHist[d];
			}

			LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);
			LIBCLOCK_PRINT(this->bIsPrintingDecodeBinTiming);
		}
		for(size_t i = 0; i < vTempEntropyField.size(); i++)
		{
			vTempEntropyField[i] = -vTempEntropyField[i] / vSum[i] + log(vSum[i]);
			// MOD-BY-LEETEN 2013/12/14:	vTempEntropyField[i] /= (ST)M_LN2;
			vTempEntropyField[i] /= (ST)log((double)uNrOfAggregatedBins);
			// MOD-BY-LEETEN 2013/12/14-END
		}
		LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);	// ADD-BY-LEETEN 01/02/2013

		LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
		_ClampToDataSize(vTempEntropyField, vEntropyField);
		_ClampBorder(vEntropyField, viLeft, viRight);
		LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);
		LIBCLOCK_PRINT(this->bIsPrintingDecodeBinTiming);
	}
	// ADD-BY-LEETEN 2013/12/01-END

	// ADD-BY-LEETEN 2013/12/07-BEGIN
	void
	_ComputeMedian(
		const vector<int>& viLeft,
		const vector<int>& viRight,
		size_t uNrOfBins,
		vector<ST>& vMedians,
		void *_Reserved = NULL
	)
	{
		LIBCLOCK_INIT(this->bIsPrintingDecodeBinTiming, __FUNCTION__);
		LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
		const size_t uDataSize = this->uDataSize;
		vector<size_t> vuDecodedLengths;

		this->_GetDecodedSize(vuDecodedLengths);

		LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);	

		// Compute the sum for all regions so we can decide the median.
		vector<ST> vSAT;	
		vector<ST> vRegionSum;	
		vector<ST> vAggregatedBin;	
		vector<ST> vRegionCumsum;	
		vector<ST> vTempMedians;	

		_DecodeAggregatedBin(0, this->UGetNrOfBins(), vSAT);

		if( vRegionSum.size() != vSAT.size() )
			vRegionSum.resize(vSAT.size());
		if( vAggregatedBin.size() != vSAT.size() )
			vAggregatedBin.resize(vSAT.size());
		if( vRegionCumsum.size() != vSAT.size() )
			vRegionCumsum.resize(vSAT.size());
		if( vTempMedians.size() != vSAT.size() )
			vTempMedians.resize(vSAT.size());
		
		vRegionCumsum.assign(vRegionCumsum.size(), (ST)0);
		vTempMedians.assign(vTempMedians.size(), (ST)0);

		// MOD-BY-LEETEN 2013/12/14:	_ComputeRegionSum(viLeft, viRight, vSAT, vRegionSum);
		_ComputeRegionSum(
			viLeft, 
			viRight, 
			vSAT, 
			vuCoefLengths,
			vRegionSum);
		// MOD-BY-LEETEN 2013/12/14-END

		/////////////// compute the SAT
		LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);	

		size_t uNrOfAggregatedBins = ( 0 == uNrOfBins ) ? this->UGetNrOfBins() : uNrOfBins ;
		size_t uBinInterval = this->UGetNrOfBins() / uNrOfAggregatedBins;
		for(size_t b = 0; b < uNrOfAggregatedBins; b++)
		{
			LIBCLOCK_INIT(this->bIsPrintingDecodeBinTiming, __FUNCTION__);
			LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
			_DecodeAggregatedBin((BT)b * uBinInterval, (BT)(b + 1) * uBinInterval - 1, vSAT);
			LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);

			// compute the local sum
			vAggregatedBin.assign(vAggregatedBin.size(), (ST)0);

			LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
			// MOD-BY-LEETEN 2013/12/14-FROM:			_ComputeRegionSum(viLeft, viRight, vSAT, vAggregatedBin);
			_ComputeRegionSum(
				viLeft, 
				viRight, 
				vSAT, 
				vuCoefLengths,
				vAggregatedBin);
			// MOD-BY-LEETEN 2013/12/14-END
			LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);

			LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
			for(size_t d = 0; d < vAggregatedBin.size(); d++)
			{
				ST dProb = vAggregatedBin[d] / vRegionSum[d];
				ST dCumProb = vRegionCumsum[d];
				if( dCumProb <= 0.5 && dCumProb + dProb > 0.5 ) {
					vTempMedians[d] = (ST)uBinInterval * ((ST)b +  (0.5 - dCumProb )/dProb); 
				}
				vRegionCumsum[d] = dCumProb + dProb;
			}

			LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);
			LIBCLOCK_PRINT(this->bIsPrintingDecodeBinTiming);
		}
		LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);

		LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
		_ClampToDataSize(vTempMedians, vMedians);
		_ClampBorder(vMedians, viLeft, viRight);
		LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);

		LIBCLOCK_PRINT(this->bIsPrintingDecodeBinTiming);
	}
	// ADD-BY-LEETEN 2013/12/07-END
};


/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
