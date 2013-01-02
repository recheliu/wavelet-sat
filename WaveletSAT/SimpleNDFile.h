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
// MOD-BY-LEETEN 01/02/2013-FROM:	_ScanHistogramForEntropy
ScanHistogramForEntropy
// MOD-BY-LEETEN 01/02/2013-END
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
	class DT,	//!< Type of the data items
	class WT	//!< Type of the wavelet coefficients.
>
class CSimpleNDFile:
	virtual public WaveletSAT::CSATSepDWTDecoder<DT, WT>
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

	// ADD-BY-LEETEN 12/29/2012-BEGIN
	void
	_DecodeSAT(
		void *_Reserved = NULL
	)
	{
		fill(this->vusCachedNextOffsets.begin(), this->vusCachedNextOffsets.end(), (unsigned short)0); 
		// fill(vusCachedBins.begin(),	vusCachedBins.end(), (unsigned short)0);
		// fill(vCachedValues.begin(),	vCachedValues.end(), (WT)0);
		valarray<DT> vSAT;
		for(size_t b = 0; b < this->UGetNrOfBins(); b++)
		{
			// LOG_VAR(b);
			_DecodeBin(b, vSAT);
		}
	}
	// ADD-BY-LEETEN 12/29/2012-END

	// ADD-BY-LEETEN 12/30/2012-BEGIN
	void
	_ComputeEntropy(
		vector<int> viLeft,
		vector<int> viRight,
		// ADD-BY-LEETEN 12/30/2012-BEGIN
		valarray<DT>& vEntropyField,
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
		// DEL-BY-LEETEN 01/02/2013:	const size_t uNrOfCoefs = this->uNrOfCoefs;
		const vector<size_t>& vuDimLengths = this->vuDimLengths;
		// DEL-BY-LEETEN 01/02/2013:	const vector<size_t>& vuCoefLengths = this->vuCoefLengths;
		// ADD-BY-LEETEN 12/30/2012-END
		////////////////////////////////////////////////////////
		// decide the offsets
		size_t uNrOfCorners = (size_t)1<<this->UGetNrOfDims();
		vector< size_t > vuCenter;	vuCenter.resize(this->UGetNrOfDims());
		#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
		for(size_t d = 0; d < this->UGetNrOfDims(); d++)
			vuCenter[d] = vuCoefLengths[d]/2;
		size_t uCenter = WaveletSAT::UConvertSubToIndex(vuCenter, vuCoefLengths);
		#else	// MOD-BY-LEETEN 01/02/2013-TO:
		for(size_t d = 0; d < this->UGetNrOfDims(); d++)
			vuCenter[d] = vuDimLengths[d]/2;
		size_t uCenter = WaveletSAT::UConvertSubToIndex(vuCenter, vuDimLengths);
		#endif	// MOD-BY-LEETEN 01/02/2013-END

		vector< long long > vllOffsets;	vllOffsets.resize(uNrOfCorners);
		vector< int > viSigns;			viSigns.resize(uNrOfCorners);
		for(size_t i = 0; i < uNrOfCorners; i++)
		{
			vector<size_t> vuSub;	vuSub.resize(this->UGetNrOfDims());
			size_t uSign = 0; 
			for(size_t 
				d = 0, j = i; 
				d < this->UGetNrOfDims(); 
				d++, j /= 2)
			{
				size_t uOffset = (0 == j % 2)?0:1;
				uSign += uOffset;
				int iSub = (int)vuCenter[d] + ((!uOffset)?viLeft[d]:viRight[d]);
				vuSub[d] = (size_t)iSub;
			}
			viSigns[i] = (0 == uSign % 2)?(-1):(+1);
			// MOD-BY-LEETEN 01/02/2013-FROM:	vllOffsets[i] = (long long)WaveletSAT::UConvertSubToIndex(vuSub, vuCoefLengths) - (long long)uCenter;
			vllOffsets[i] = (long long)WaveletSAT::UConvertSubToIndex(vuSub, vuDimLengths) - (long long)uCenter;
			// MOD-BY-LEETEN 01/02/2013-END
		}
		LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);	// ADD-BY-LEETEN 01/02/2013

		/////////////// compute the SAT
		LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);	// ADD-BY-LEETEN 01/02/2013
		valarray<DT> vSAT;
		#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
		valarray<DT> vLocalHist;	vLocalHist.resize(uNrOfCoefs);
		valarray<DT> vTempEntropyField;	vTempEntropyField.resize(uNrOfCoefs);
		valarray<DT> vSum;			vSum.resize(uNrOfCoefs);
		#else	// MOD-BY-LEETEN 01/02/2013-TO:
		valarray<DT> vLocalHist;	
		valarray<DT> vTempEntropyField;	
		valarray<DT> vSum;			
		#endif	// MOD-BY-LEETEN 01/02/2013-END
		for(size_t b = 0; b < this->UGetNrOfBins(); b++)
		{
			// ADD-BY-LEETEN 01/02/2013-BEGIN
			LIBCLOCK_INIT(this->bIsPrintingDecodeBinTiming, __FUNCTION__);
			LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
			// ADD-BY-LEETEN 01/02/2013-END

			_DecodeBin((unsigned short)b, vSAT);

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
			vLocalHist = (DT)0;
			LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);	// ADD-BY-LEETEN 01/02/2013
			#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
			for(size_t i = 0; i < uNrOfCorners; i++)
				vLocalHist += vSAT.shift((int)vllOffsets[i]) * (DT)viSigns[i];
			// compute the entropy
			vLocalHist = vLocalHist.apply(_ClampNegative);
			vTempEntropyField += vLocalHist.apply(_ScanHistogramForEntropy);
			vSum += vLocalHist;
			#else	// MOD-BY-LEETEN 01/02/2013-TO:
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
					vLocalHist[d] += vSAT[d + vllOffsets[i]] * (DT)viSigns[i];
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
			#endif	// MOD-BY-LEETEN 01/02/2013-END

			// ADD-BY-LEETEN 01/02/2013-BEGIN
			LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);
			LIBCLOCK_PRINT(this->bIsPrintingDecodeBinTiming);
			// ADD-BY-LEETEN 01/02/2013-END
		}
		vTempEntropyField = -vTempEntropyField / vSum + log(vSum);
		vTempEntropyField /= (DT)log(2.0);
		LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);	// ADD-BY-LEETEN 01/02/2013

		#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
		// only keep the entropy field within the data range
		if( uDataSize != vEntropyField.size() )
			vEntropyField.resize(uDataSize);
		vector<size_t> vuSub;
		for(size_t d = 0; d < uDataSize; d++)
		{
			vector<size_t> vuSub;
			WaveletSAT::_ConvertIndexToSub(d, vuSub, vuDimLengths);
			// ADD-BY-LEETEN 12/30/2012-BEGIN
			bool bIsNearBorder = false;
			for(size_t dim = 0; dim < this->UGetNrOfDims(); dim++)
				if( 0 > (int)vuSub[dim] + viLeft[dim] || 
						(int)vuSub[dim] + viLeft[dim] >= vuDimLengths[dim] ||
					0 > (int)vuSub[dim] + viRight[dim] || 
						(int)vuSub[dim] + viRight[dim] >= vuDimLengths[dim] )
				{
					bIsNearBorder = true;
					break;
				}

			if( bIsNearBorder )
			{
				vEntropyField[d] = (DT)0;
				continue;
			}
			// ADD-BY-LEETEN 12/30/2012-END
			vEntropyField[d] = vTempEntropyField[WaveletSAT::UConvertSubToIndex(vuSub, vuCoefLengths)];
		}	
		#else	// MOD-BY-LEETEN 01/02/2013-TO:
		LIBCLOCK_BEGIN(this->bIsPrintingDecodeBinTiming);
		_ClampToDataSize(vTempEntropyField, vEntropyField);
		_ClampBorder(vEntropyField, viLeft, viRight);
		LIBCLOCK_END(this->bIsPrintingDecodeBinTiming);
		LIBCLOCK_PRINT(this->bIsPrintingDecodeBinTiming);
		#endif	// MOD-BY-LEETEN 01/02/2013-END
	}
	// ADD-BY-LEETEN 12/30/2012-END
};

