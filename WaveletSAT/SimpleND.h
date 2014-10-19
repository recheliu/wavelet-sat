#pragma once

#include <exception>
#include <unordered_map>
using namespace std;
#include "contourspectrum.h"

#include "WaveletSATEncoder.h"
#if	WITH_CUDA
#include "WaveletSATGPUEncoder.h"	
#endif	// #if	WITH_CUDA
#include "SATFileEncoder.h"
#if WITH_NETCDF 
#include <netcdf.h>
#endif

template<
	typename DT,							//!< Type of the data
	typename ST = WaveletSAT::typeSum,		//!< Type of the sum
	typename BT = WaveletSAT::typeBin,		//!< Type of the bin
	typename WT = WaveletSAT::typeWavelet	//!< Type of the wavelet coefficientsd
>
class CSimpleND:
	#if	!WITH_SAT_FILE	
	#if	!WITH_CUDA		
	virtual public WaveletSAT::CWaveletSATEncoder<DT, ST, BT, WT>

	#else	// #if	!WITH_CUDA	
	virtual public WaveletSAT::CWaveletSATGPUEncoder<DT, ST, BT, WT>
	#endif	// #if	!WITH_CUDA	

	#else	// #if	!WITH_SAT_FILE	
	virtual public WaveletSAT::CSATFileEncoder<DT, ST, BT>
	#endif	// #if	!WITH_SAT_FILE	
{
	DT valueMin, valueMax;
	bool bIsWithContourSpectrum;
	double dMinSum;
	const vector<DT>* pvData;

	vector<double> vdBinEdges;
public:
	enum EParameter
	{
		PARAMETER_BEGIN = 0x0F00,

		WITH_CONTOUR_SPECTRUM,

		PARAMETER_END
	};

	virtual
	void
	_SetData(
		const vector<DT>* pvData,
		void* _Reserved = NULL
	)
	{
		this->pvData = pvData;
	}

	virtual
	void
	_SetInteger(
		int eName,
		long lValue,
		void* _Reserved = NULL
	)
	{
		switch(eName)
		{
		case WITH_CONTOUR_SPECTRUM:
			bIsWithContourSpectrum = (lValue)?true:false;
			break;
		}

		#if	!WITH_SAT_FILE	
		#if	!WITH_CUDA		
		WaveletSAT::CWaveletSATEncoder<DT, ST, BT, WT>::_SetInteger(eName, lValue);
		#else	// #if	!WITH_CUDA	
		WaveletSAT::CWaveletSATGPUEncoder<DT, ST, BT, WT>::_SetInteger(eName, lValue);
		#endif	// #if	!WITH_CUDA	
		#else	// #if	!WITH_SAT_FILE	
		WaveletSAT::CSATFileEncoder<DT, ST, BT>::_SetInteger(eName, lValue);
		#endif	// #if	!WITH_SAT_FILE	
	}

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
	  const size_t& uNrOfBins = this->uNrOfBins; 
		if( this->bIsWithContourSpectrum && 
			(3 == vuPos.size() || 2 == vuPos.size() ) )
		{
			if(vdBinEdges.empty())
			{
				vdBinEdges.resize(uNrOfBins+1);
				double dBinInterval = (double)(valueMax - valueMin)/(double)uNrOfBins;
				for(size_t b = 0; b < vdBinEdges.size(); b++)
					vdBinEdges[b] = dBinInterval * (double)b + (double)valueMin;
			}

			size_t uOrientation = 0;
			bool bIsBorder = false;
			for(size_t d = 0; d < vuPos.size(); d++)
			{
				if( vuPos[d] % 2 ) 
					uOrientation = 1 - uOrientation;

				if( !vuPos[d] ) 
				{
					bIsBorder = true;
					break;
				}
			}

			vpBins.clear();
			if( !bIsBorder )
			{
				vector<pair<double, glm::dvec4> > vCorners;
				vector<size_t> vuCorner;
				vuCorner.resize(vuPos.size());
				vector<size_t> vuOffset, vuOffsetMax;
				vuOffsetMax.assign(vuPos.size(), 2);
				size_t uNrOfCorners = 1 << vuPos.size();
				for(size_t c = 0; c < uNrOfCorners; c++)
				{
					WaveletSAT::_ConvertIndexToSub(c, vuOffset, vuOffsetMax);
					glm::dvec4 vd4;
					for(size_t d = 0; d < vuPos.size(); d++)
					{
						vuCorner[d] = vuPos[d] - (1 - vuOffset[d]);
						vd4[d] = (double)vuCorner[d];
					}

					for(size_t d = vuPos.size(); d < 4; d++)
						vd4[d] = 0.0;

					double dScalar = (double)(*pvData)[WaveletSAT::UConvertSubToIndex(vuCorner, this->vuDimLengths)];
					vCorners.push_back(make_pair<double, glm::dvec4>(dScalar, vd4));
				}
				vector<double> vdSpectrum;
				vdSpectrum.assign(uNrOfBins, 0.0);
				size_t uFirstBin = uNrOfBins;
				switch(vuPos.size())
				{
				case 2:
					ContourSpectrum::_ComputeFor2DCell
					(
						vdBinEdges,
						vCorners,
						uFirstBin,
						vdSpectrum,
						uOrientation
					);
					break;

				case 3:
					ContourSpectrum::_ComputeFor3DCell
					(
						vdBinEdges,
						vCorners,
						uFirstBin,
						vdSpectrum,
						uOrientation
					);
					break;
				}
				for(size_t b = uFirstBin; b < uNrOfBins; b++)
				{
					if( !vdSpectrum[b] )
						break;
					vpBins.push_back(pair<BT, ST>((BT)b, (ST)vdSpectrum[b]));
					dMinSum = min(dMinSum, vdSpectrum[b]);
				}
			}
		}
		else
		{
			DT clampedValue = min(max(value, valueMin), valueMax);
			size_t uBin = (size_t)floorf((float)(uNrOfBins * (clampedValue - valueMin))/(float)(valueMax - valueMin));
			uBin = min(uBin, uNrOfBins - 1);
			vpBins.clear();
			vpBins.push_back(pair<BT, ST>((BT)uBin, (ST)1));
		}	
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

	#if !WITH_SAT_FILE
	double 
	DGetThreshold
	(
		void *_Reserved = NULL
	)
	{
		double dWeight = 1.0;
		if( bIsWithContourSpectrum )
			if(	dMinSum < HUGE_VAL )
				dWeight = dMinSum;
			else
				throw std::range_error("Invalide dMinSum.");

		return dWeight * dWaveletThreshold;
	}
	#endif	// #if	!WITH_SAT_FILE

	CSimpleND<DT, ST, BT, WT>()
	{
		bIsWithContourSpectrum = false;
		dMinSum = HUGE_VAL;
		pvData = NULL;
	}
};


/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
