#pragma once

#include	<algorithm>	// ADD-BY-LEETEN 03/17/2013

#include <unordered_map>	
#include <vector>
using namespace std;
#include <math.h>

#include "HeaderBase.h"	
#include "DWT.h"
#include "liblog.h"	

namespace WaveletSAT
{
	//! SepDWT coefficients per basis of all bins
	template<
		typename WT = typeWavelet,	//!< Type of the wavelet coefficientsd
		typename BT = typeBin		//!< Type of the bin
	>
	class CSepDWTPoolBase
		:virtual public CBase	// ADD-BY-LEETEN 12/30/2012
	{
	protected:	
		// the pool of coefficients 
		//! Size of this coef. block. 
		vector<size_t> vuLengths;

		//! Size of this coef. block, which is equal to the product of vuLengths
		size_t uSize;

		//! A flag to indicate whether this coef. block is sparse.
		bool bIsSparse;

		vector<size_t> vuDataDimLengths;

		vector<size_t> vuWaveletLengths;

		WT	 dWaveletWeight;

		size_t uNrOfBins;

		typedef pair<size_t, vector<WT>*> CFullArray;
		typedef unordered_map<size_t,  CFullArray*> CFullArrays;
		CFullArrays* pcFullArrays;
		
		//! Max # that each coefficient is updated
		size_t uMaxCount;

	protected:
		void
		_GetAtFullArray
		(
			const BT& Bin,
			const size_t& uIndex,
			WT& Value,
			void* _Reserved = NULL
		) 
		{
			CFullArrays::iterator iterFullArrays = this->pcFullArrays->find(uIndex);
			if( pcFullArrays->end() 
					 == iterFullArrays || 
				NULL == iterFullArrays->second ||
				NULL == iterFullArrays->second->second 
				)
			{
				Value = WT(0);
				return;
			}
			vector<WT>& vFullArray = *iterFullArrays->second->second;
			Value = vFullArray[Bin];
		}

		size_t UGetNonOccupiedVolume(
			const size_t& uIndex, 
			void *_Reserved = NULL
			)
		{
			static vector<size_t> vuSubs;
			_ConvertIndexToSub(uIndex, vuSubs, vuLengths);

			size_t uOccupiedVolume = 1;
			for(size_t d = 0; d < vuSubs.size(); d++) 
			{
				size_t uMaxOccupiedLength = (vuSubs[d] + 1) * vuWaveletLengths[d];
				size_t uOccupiedLength = vuWaveletLengths[d];
				if( uMaxOccupiedLength > this->vuDataDimLengths[d] )
				{
					size_t uMinOccupiedLength = vuSubs[d] * vuWaveletLengths[d];
					if( uMinOccupiedLength >= this->vuDataDimLengths[d] )
					{
						uOccupiedVolume = 0;
						break;
					}
					else
						uOccupiedLength = this->vuDataDimLengths[d] - uMinOccupiedLength;
				}
				uOccupiedVolume *= uOccupiedLength;
			}
			return uMaxCount - uOccupiedVolume;
		}

	public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		bool 
		BIsSparse
		(
			void* _Reserved = NULL
		)
		{
			return bIsSparse;
		}

		void
		_SetDataDimLengths
		(
			const vector<size_t>& vuDataDimLengths,
			void* _Reserved = NULL
			)
		{
			this->vuDataDimLengths.assign(vuDataDimLengths.begin(), vuDataDimLengths.end());
		}

		void
		_SetWaveletLengths
		(
			const vector<size_t>& vuWaveletLengths,
			void* _Reserved = NULL
			)
		{
			this->vuWaveletLengths.assign(vuWaveletLengths.begin(), vuWaveletLengths.end());
		}

		void
		_SetWaveletWeight(
			const WT& dWaveletWeight,
			void* _Reserved = NULL
			)
		{
			this->dWaveletWeight = dWaveletWeight;	
		}

		void
		_Set
		(
			const BT& uNrOfBins,
			const vector<size_t>& vuLengths,
			size_t uMaxCount,	
			bool bIsSparse,
			void* _Reserved = NULL
		)
		{
			this->vuLengths.clear();
			this->vuLengths.resize(vuLengths.size());
			this->uSize = 1;
			for(size_t d = 0; d < vuLengths.size(); d++)
			{
				this->vuLengths[d] = vuLengths[d];
				this->uSize *= vuLengths[d];
			}
			this->bIsSparse = bIsSparse;
			this->uNrOfBins = uNrOfBins;	
			if( !bIsSparse )
			{
				this->pcFullArrays = new CFullArrays();
			}
			this->uMaxCount = uMaxCount;	
		}

		void
		_IncreaseCount
		(
			size_t uExtraCount,
			void *_Reserved = NULL
		)
		{
			uMaxCount += uExtraCount;
		}

		void
		_AddAtFullArray
		(
			const BT& usBin,
			const size_t& uIndex,
			const WT& Value,
			const size_t& uCount,	
			void* _Reserved = NULL
		) 
		{
			CFullArrays::iterator iterFullArrays = this->pcFullArrays->find(uIndex);
			CFullArray* pcFullArray = NULL;
			if( iterFullArrays != pcFullArrays->end() )
			{
				pcFullArray = iterFullArrays->second;
			}
			else
			{
				pcFullArray = new CFullArray();
				this->pcFullArrays->insert(pair<size_t, CFullArray*>(uIndex, pcFullArray));
			}

			if( !pcFullArray->second ) 
			{
				pcFullArray->first = UGetNonOccupiedVolume(uIndex);

				pcFullArray->second = new vector<WT>();
				pcFullArray->second->assign(uNrOfBins, (WT)0);
			} 
			pcFullArray->first += uCount;
			vector<WT>& vec = *pcFullArray->second;
			vec[usBin] += Value;

			if( !uCount ) 
				return;
		}

		void
		_Finalize
		(
			WT WaveletWeight,
			void* _Reserved = NULL
		)
		{
		}


		CSepDWTPoolBase()
		{
			this->pcFullArrays = NULL;
		}

		virtual	
		~CSepDWTPoolBase()
		{
			if( this->pcFullArrays ) {
				for(CFullArrays::iterator 
						iterFullArrays = pcFullArrays->begin(); 
					iterFullArrays != pcFullArrays->end(); 
					iterFullArrays++) 
				{
					if( iterFullArrays->second ) 
					{
						if( iterFullArrays->second->second ) 
							delete iterFullArrays->second->second;
						delete iterFullArrays->second;
					}
				}
			}
		}
	};
}
		
