#pragma once

#include	<algorithm>	// ADD-BY-LEETEN 03/17/2013

#include <unordered_map>	
#include <vector>
using namespace std;
#include <math.h>

#include "HeaderBase.h"	
#include "DWT.h"
#include "SATSepDWTNetCDF.h"	
#include "SepDWTPoolBase.h"
#include "liblog.h"	

namespace WaveletSAT
{
	//! SepDWT coefficients per basis of all bins
	template<
		typename WT = typeWavelet,	//!< Type of the wavelet coefficientsd
		typename BT = typeBin		//!< Type of the bin
	>
	class CSepDWTDecoderPool
		:virtual public CSepDWTPoolBase<WT, BT>	// ADD-BY-LEETEN 12/30/2012
	{
	protected:	
		typedef vector<pair<BT, WT>> CDecodingSparseArray;
		typedef unordered_map< size_t, CDecodingSparseArray* > CDecodingSparseArrays;
		CDecodingSparseArrays *pcDecodingSparseArrays;

	protected:
		void
		_AppendDecodingSparseArray
		(
			const BT& usBin,
			const size_t& uIndex,
			const WT& Value,
			const size_t& uCount,
			void* _Reserved = NULL
		) 
		{
			CDecodingSparseArrays::iterator iterSparseArrays = this->pcDecodingSparseArrays->find(uIndex);
			CDecodingSparseArray* pcSparseArray = NULL;
			if( iterSparseArrays != pcDecodingSparseArrays->end() )
			{
				pcSparseArray = iterSparseArrays->second;
			}
			else
			{
				pcSparseArray = new CDecodingSparseArray();
				this->pcDecodingSparseArrays->insert(pair<size_t, CDecodingSparseArray*>(uIndex, pcSparseArray));
			}

			pcSparseArray->push_back(pair<BT, WT>(usBin, Value));
		}

		void
		_GetAtDecodingSparseArray
		(
			const BT& usOffset,
			const size_t& uIndex,
			BT& Bin,
			WT& Value,
			void* _Reserved = NULL
		) 
		{
			CDecodingSparseArrays::iterator iterDecodingSparseArrays = pcDecodingSparseArrays->find(uIndex);
			if( pcDecodingSparseArrays->end()
					 != iterDecodingSparseArrays &&
				NULL != iterDecodingSparseArrays->second )
			{
				vector<pair<BT, WT>>& vpair = *iterDecodingSparseArrays->second;
				if( (size_t)usOffset < vpair.size() )
				{
					Bin = vpair[usOffset].first;
					Value = vpair[usOffset].second;
				}
			}
		}

	public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/
		void
		_Copy(
			size_t uIndex,
			size_t uNrOfBins,
			const CSATSepDWTNetCDF::TYPE_COEF_BIN	*pBins,
			const CSATSepDWTNetCDF::TYPE_COEF_VALUE	*pValues,
			void *_Reserved = NULL
		)
		{
			for(size_t b = 0; b < uNrOfBins; b++)
			{
				CSATSepDWTNetCDF::TYPE_COEF_VALUE	Value = pValues[b];
				CSATSepDWTNetCDF::TYPE_COEF_BIN		Bin =	pBins[b];
				if( Value )
				{
					if( !bIsSparse )
						this->_AddAtFullArray(Bin, uIndex, Value, 0);
					else
						this->_AppendDecodingSparseArray(Bin, uIndex, Value, 0);
				}	
			}
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
			CSepDWTPoolBase::_Set(
				uNrOfBins, 
				vuLengths,
				uMaxCount,
				bIsSparse,
				_Reserved
			);

			if( bIsSparse )
			{	
				this->pcDecodingSparseArrays = new CDecodingSparseArrays();
			}
		}

		void
		_GetAtOffset
		(
			const BT& usOffset,
			const vector<size_t>& vuSubs,
			size_t& uIndex,
			BT& Bin,
			WT& Value,
			void* _Reserved = NULL
		)
		{
			if( vuSubs.size() )
				uIndex = UConvertSubToIndex(vuSubs, vuLengths);

			if( !bIsSparse )
			{
				// full
				Bin = usOffset;
				_GetAtFullArray(Bin, uIndex, Value);
			}
			else
			{
				_GetAtDecodingSparseArray(usOffset, uIndex, Bin, Value);
			}
		}

		void
		_GetCoefSparse
		(
			size_t uIndex,
			vector< pair<BT, WT> >& vpairCoefs,
			void* _Reserved = NULL
		) const
		{
			vpairCoefs.clear();
			if( !bIsSparse )
			{
				CFullArrays::iterator iterFullArrays = pcFullArrays->find(uIndex);
				if( pcFullArrays->end() 
						 != iterFullArrays &&
					NULL != iterFullArrays->second &&
					NULL != iterFullArrays->second->second
					) 
				{
					vector<WT>& vFullArray = *iterFullArrays->second->second;
					for(size_t b = 0; b < uNrOfBins ; b++)
					{
						WT Coef = vFullArray[b];
						if( Coef )
							vpairCoefs.push_back(pair<BT, WT>((BT)b, Coef));
					}
				}
			}
			else
			{
				CDecodingSparseArrays::iterator iterSparseArrays = pcDecodingSparseArrays->find(uIndex);
				if( pcDecodingSparseArrays->end() 
						 != iterSparseArrays && 
					NULL != iterSparseArrays->second )
				{
					const vector< pair<BT, WT> >& vpairSparse = *iterSparseArrays->second;
					vpairCoefs.resize(vpairSparse.size());
					copy(vpairSparse.begin(), vpairSparse.end(), vpairCoefs.begin());
				}
			}	
		}
		
		virtual
		void
		_GetCoefSparse
		(
			const vector<size_t> vuSubs,
			vector< pair<BT, WT> >& vpairCoefs,
			void* _Reserved = NULL
		) const
		{
			size_t uIndex = UConvertSubToIndex(vuSubs, vuLengths);
			_GetCoefSparse
			(
				uIndex,
				vpairCoefs
				);
		}

		CSepDWTDecoderPool():
			CSepDWTPoolBase()
		{
		}

		// ADD-BY-LEETEN 11/12/2012-BEGIN
		virtual
		~CSepDWTDecoderPool()
		{
			if( this->pcDecodingSparseArrays ) {
				for(CDecodingSparseArrays::iterator 
						iterSparseArrays = pcDecodingSparseArrays->begin(); 
					iterSparseArrays != pcDecodingSparseArrays->end(); 
					iterSparseArrays++) 
				{
					if( iterSparseArrays->second ) 
					{
						delete iterSparseArrays->second;
					}
				}
			}
		}
		// ADD-BY-LEETEN 11/12/2012-END
	};
}
		
