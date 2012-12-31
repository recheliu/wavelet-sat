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

// ADD-BY-LEETEN 12/30/2012-BEGIN
//! Scan the histogra bin to update the entropy
template<class DT>	//!< Type of the data items
DT
_ScanHistogramForEntropy
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
// MOD-BY-LEETEN 12/29/2012-FROM:	template<class DT>
template<
	class DT,	//!< Type of the data items
	class WT	//!< Type of the wavelet coefficients.
>
// MOD-BY-LEETEN 12/29/2012-END
class CSimpleNDFile:
	// MOD-BY-LEETEN 12/25/2012-FROM:	public WaveletSAT::CSATSepDWTFile<DT>
	// MOD-BY-LEETEN 12/29/2012-FROM:	public WaveletSAT::CSATSepDWTOutOfCore<DT>
	// MOD-BY-LEETEN 12/30/2012-FROM:	public WaveletSAT::CSATSepDWTOutOfCore<DT, WT>
	virtual public WaveletSAT::CSATSepDWTOutOfCore<DT, WT>
	// MOD-BY-LEETEN 12/30/2012-END
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

	// ADD-BY-LEETEN 12/29/2012-BEGIN
	void
	_DecodeSAT(
		void *_Reserved = NULL
	)
	{
	  #if 0 	  // MOD-BY-LEETEN 12/30/2012-FROM:
		fill(vusCachedNextOffsets.begin(), vusCachedNextOffsets.end(), (unsigned short)0);
		#else // MOD-BY-LEETEN 12/30/2012-TO:
		fill(this->vusCachedNextOffsets.begin(), this->vusCachedNextOffsets.end(), (unsigned short)0); 
		#endif // MOD-BY-LEETEN 12/30/2012-END
		// fill(vusCachedBins.begin(),	vusCachedBins.end(), (unsigned short)0);
		// fill(vCachedValues.begin(),	vCachedValues.end(), (WT)0);
		// MOD-BY-LEETEN 12/30/2012-FROM:		vector<DT> vSAT;
		valarray<DT> vSAT;
		// MOD-BY-LEETEN 12/30/2012-END
		#if 0 		// MOD-BY-LEETEN 12/30/2012-FROM:
		for(size_t b = 0; b < UGetNrOfBins(); b++)
		  #else // MOD-BY-LEETEN 12/30/2012-TO:
		for(size_t b = 0; b < this->UGetNrOfBins(); b++)
		  #endif // MOD-BY-LEETEN 12/30/2012-END
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
	  // ADD-BY-LEETEN 12/30/2012-BEGIN
		const size_t uDataSize = this->uDataSize;
		const size_t uNrOfCoefs = this->uNrOfCoefs;
		const vector<size_t>& vuDimLengths = this->vuDimLengths;
		const vector<size_t>& vuCoefLengths = this->vuCoefLengths;
		// ADD-BY-LEETEN 12/30/2012-END
		////////////////////////////////////////////////////////
		// decide the offsets
		#if 0 		// MOD-BY-LEETEN 12/30/2012-FROM:
		size_t uNrOfCorners = 1<<UGetNrOfDims();
		vector< size_t > vuCenter;	vuCenter.resize(UGetNrOfDims());
		for(size_t d = 0; d < UGetNrOfDims(); d++)
			vuCenter[d] = vuCoefLengths[d]/2;
		#else // MOD-BY-LEETEN 12/30/2012-TO:
		size_t uNrOfCorners = (size_t)1<<this->UGetNrOfDims();
		vector< size_t > vuCenter;	vuCenter.resize(this->UGetNrOfDims());
		for(size_t d = 0; d < this->UGetNrOfDims(); d++)
			vuCenter[d] = vuCoefLengths[d]/2;
		#endif // MOD-BY-LEETEN 12/30/2012-END
		size_t uCenter = WaveletSAT::UConvertSubToIndex(vuCenter, vuCoefLengths);

		vector< long long > vllOffsets;	vllOffsets.resize(uNrOfCorners);
		vector< int > viSigns;			viSigns.resize(uNrOfCorners);
		for(size_t i = 0; i < uNrOfCorners; i++)
		{
		  #if 0 		  // MOD-BY-LEETEN 12/30/2012-FROM:
			vector<size_t> vuSub;	vuSub.resize(UGetNrOfDims());
			#else // MOD-BY-LEETEN 12/30/2012-TO:
			vector<size_t> vuSub;	vuSub.resize(this->UGetNrOfDims());
			#endif // MOD-BY-LEETEN 12/30/2012-END
			size_t uSign = 0; 
			for(size_t 
				d = 0, j = i; 
			    #if 0 			    // MOD-BY-LEETEN 12/30/2012-FROM:
				d < UGetNrOfDims(); 
			    #else // MOD-BY-LEETEN 12/30/2012-TO:
				d < this->UGetNrOfDims(); 
			    #endif // MOD-BY-LEETEN 12/30/2012-END
				d++, j /= 2)
			{
				size_t uOffset = (0 == j % 2)?0:1;
				uSign += uOffset;
				int iSub = (int)vuCenter[d] + ((!uOffset)?viLeft[d]:viRight[d]);
				vuSub[d] = (size_t)iSub;
			}
			viSigns[i] = (0 == uSign % 2)?(-1):(+1);
			vllOffsets[i] = (long long)WaveletSAT::UConvertSubToIndex(vuSub, vuCoefLengths) - (long long)uCenter;
		}

		/////////////// compute the SAT
		valarray<DT> vSAT;
		valarray<DT> vLocalHist;	vLocalHist.resize(uNrOfCoefs);
		valarray<DT> vTempEntropyField;	vTempEntropyField.resize(uNrOfCoefs);
		valarray<DT> vSum;			vSum.resize(uNrOfCoefs);
		#if 0 		// MOD-BY-LEETEN 12/30/2012-FROM:
		for(size_t b = 0; b < UGetNrOfBins(); b++)
		  #else // MOD-BY-LEETEN 12/30/2012-TO:
		for(size_t b = 0; b < this->UGetNrOfBins(); b++)
		  #endif // MOD-BY-LEETEN 12/30/2012-END
		{
			// MOD-BY-LEETEN 12/31/2012-FROM:			_DecodeBin(b, vSAT);
			_DecodeBin((unsigned short)b, vSAT);
			// MOD-BY-LEETEN 12/31/2012-END

			// compute the local sum
			vLocalHist = (DT)0;
			for(size_t i = 0; i < uNrOfCorners; i++)
				// MOD-BY-LEETEN 12/31/2012-FROM:	vLocalHist += vSAT.shift(vllOffsets[i]) * (DT)viSigns[i];
				vLocalHist += vSAT.shift((int)vllOffsets[i]) * (DT)viSigns[i];
				// MOD-BY-LEETEN 12/31/2012-END

			// compute the entropy
			vLocalHist = vLocalHist.apply(_ClampNegative);
			vTempEntropyField += vLocalHist.apply(_ScanHistogramForEntropy);
			vSum += vLocalHist;
		}
		vTempEntropyField = -vTempEntropyField / vSum + log(vSum);
		vTempEntropyField /= (DT)log(2.0);

		// only keep the entropy field within the data range
		// MOD-BY-LEETEN 12/30/2012-FROM:		valarray<DT> vEntropyField;	vEntropyField.resize(uDataSize);
		if( uDataSize != vEntropyField.size() )
			vEntropyField.resize(uDataSize);
		// MOD-BY-LEETEN 12/30/2012-END
		vector<size_t> vuSub;
		for(size_t d = 0; d < uDataSize; d++)
		{
			vector<size_t> vuSub;
			WaveletSAT::_ConvertIndexToSub(d, vuSub, vuDimLengths);
			// ADD-BY-LEETEN 12/30/2012-BEGIN
			bool bIsNearBorder = false;
			#if 0 			// MOD-BY-LEETEN 12/30/2012-FROM:
			for(size_t dim = 0; dim < UGetNrOfDims(); dim++)
			  #else // MOD-BY-LEETEN 12/30/2012-TO:
			for(size_t dim = 0; dim < this->UGetNrOfDims(); dim++)
			  #endif // MOD-BY-LEETEN 12/30/2012-END
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
	}
	// ADD-BY-LEETEN 12/30/2012-END
};

