#pragma once

#include	<algorithm>	

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
		:virtual public CSepDWTPoolBase<WT, BT>	
	{
	protected:	
		typedef vector<pair<BT, WT>> CDecodingSparseArray;
		#if	WITHOUT_FULL_ARRAYS	
			typedef unordered_map< size_t, CDecodingSparseArray* > CDecodingSparseArrays;
			CDecodingSparseArrays *pcDecodingSparseArrays;
		#else	// #if	WITHOUT_FULL_ARRAYS	
			typedef vector< CDecodingSparseArray* > CDecodingSparseArrays;
			CDecodingSparseArrays vcDecodingSparseArrays;
		#endif	// #if	WITHOUT_FULL_ARRAYS	

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
			#if	WITHOUT_FULL_ARRAYS		
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
			#else	// #if	WITHOUT_FULL_ARRAYS	
			CDecodingSparseArray* pcSparseArray = vcDecodingSparseArrays[uIndex];
			if( !pcSparseArray )
			{
				pcSparseArray = new CDecodingSparseArray();
				vcDecodingSparseArrays[uIndex] = pcSparseArray;
			}
			#endif	// #if	WITHOUT_FULL_ARRAYS	

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
			#if	WITHOUT_FULL_ARRAYS		
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
			#else	// #if	WITHOUT_FULL_ARRAYS	
			CDecodingSparseArray* pcSparseArray = pcDecodingSparseArrays[uIndex];
			if( pcSparseArray )
			{
				vector<pair<BT, WT>>& vpair = *pcSparseArray->second;
				if( (size_t)usOffset < vpair.size() )
				{
					Bin = vpair[usOffset].first;
					Value = vpair[usOffset].second;
				}
			}
			#endif	// #if	WITHOUT_FULL_ARRAYS	
		}

	public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/
		void
		_Finalize
		(
			WT WaveletWeight,
			void* _Reserved = NULL
		)
		{
			if( bIsSparse ) {
				#if	WITHOUT_FULL_ARRAYS		
				for(CDecodingSparseArrays::iterator 
						iterSparseArrays = pcDecodingSparseArrays->begin(); 
					iterSparseArrays != pcDecodingSparseArrays->end(); 
					iterSparseArrays++) 
				{
					if( iterSparseArrays->second ) 
					{
						CDecodingSparseArray& sparseArray = *iterSparseArrays->second;
				#else	// #if	WITHOUT_FULL_ARRAYS	
				for(CDecodingSparseArrays::iterator 
						iterSparseArrays = vcDecodingSparseArrays.begin(); 
					iterSparseArrays != vcDecodingSparseArrays.end(); 
					iterSparseArrays++) 
				{
					if( *iterSparseArrays ) 
					{
						CDecodingSparseArray& sparseArray = *(*iterSparseArrays);
				#endif	// #if	WITHOUT_FULL_ARRAYS	
						sort(sparseArray.begin(), sparseArray.end());
		
						#if	!WITHOUT_BIN_AGGREGATION
						WT dCumsum = 0.0;
						for(vector<pair<BT, WT>>::iterator
								iter = sparseArray.begin();
							iter != sparseArray.end();
							iter++) 
						{
							dCumsum += iter->second;
							iter->second = dCumsum;
						}
						#endif	// #if	!WITHOUT_BIN_AGGREGATION
					}
				}
			}

			#if	!WITHOUT_BIN_AGGREGATION
			else
			{
				for(CFullArrays::iterator 
						iterFullArrays = pcFullArrays->begin(); 
					iterFullArrays != pcFullArrays->end(); 
					iterFullArrays ++) 
				{
					if( pcFullArrays->end() 
							 != iterFullArrays &&
						NULL != iterFullArrays->second &&
						NULL != iterFullArrays->second->second
						) 
					{
						vector<WT>& vdValues = *iterFullArrays->second->second;
						WT dCumsum = 0.0;
						for(vector<WT>::iterator 
								iterValues = vdValues.begin();
							iterValues != vdValues.end();
							iterValues++) 
						{
							dCumsum += *iterValues;
							*iterValues = dCumsum;
						}
					}
				}
			}
			#endif	// #if	!WITHOUT_BIN_AGGREGATION
		}

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
				#if	WITHOUT_FULL_ARRAYS		
					this->pcDecodingSparseArrays = new CDecodingSparseArrays();
				#else	// #if	WITHOUT_FULL_ARRAYS	
					vcDecodingSparseArrays.resize(UGetProduct(vuLengths));
				#endif	// #if	WITHOUT_FULL_ARRAYS	
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

		#if	!WITHOUT_BIN_AGGREGATION
		static 
		bool 
		BCompareBin
		(
			const pair<BT, WT>& pairBin1,
			const pair<BT, WT>& pairBin2
		)
		{
			return ( pairBin1.first < pairBin2.first ) ;
		}

		// Return the cumsums till the left edge of the specified bins.
		virtual	
		void
		_GetCumsums
		(
			size_t uIndex,
			const vector<BT>& vsBins,
			vector<WT>& udValues,
			void* _Reserved = NULL
		) const 
		{
			udValues.clear();
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
					for(vector<BT>::const_iterator 
							iterBins = vsBins.begin();
						iterBins != vsBins.end();
						iterBins ++ ) 
					{
						BT usBin = *iterBins;
						WT dCumsum = ( !usBin ) ? 0.0 : vFullArray[min((size_t)usBin - 1, vFullArray.size()-1)];
						udValues.push_back(dCumsum);
					}
				}
			}
			else
			{
				#if	WITHOUT_FULL_ARRAYS		
				CDecodingSparseArrays::iterator iterSparseArrays = pcDecodingSparseArrays->find(uIndex);
				if( pcDecodingSparseArrays->end() 
						 != iterSparseArrays && 
					NULL != iterSparseArrays->second )
				{
					const vector< pair<BT, WT> >& vpairSparse = *iterSparseArrays->second;
				#else	// #if	WITHOUT_FULL_ARRAYS	
				if( vcDecodingSparseArrays[uIndex] )
				{
					const vector< pair<BT, WT> >& vpairSparse = *vcDecodingSparseArrays[uIndex];
				#endif	// #if	WITHOUT_FULL_ARRAYS	

					for(vector<BT>::const_iterator 
							iterBins = vsBins.begin();
						iterBins != vsBins.end();
						iterBins ++ ) 
					{
						BT usBin = *iterBins;
						WT dCumsum = 0.0;
						/*
						// The code after the comment is equivalent to the pseudo code here.

						dCumsum = 0.0;
						for(vector< pair<BT, WT> >::const_iterator 
								iterLowerBound = vpairSparse.begin();
							iterLowerBound != vpairSparse.end();
							iterLowerBound ++) 
						{
							if( usBin <= iterLowerBound->first ) 
							{
								break;
							}
							dCumsum = iterLowerBound->second;
						}
						*/
						pair<BT, WT> pairValue = make_pair<BT, WT>(usBin, (WT)0.0);	
						vector< pair<BT, WT> >::const_iterator iterLowerBound = 
							lower_bound(vpairSparse.begin(), vpairSparse.end(), pairValue, BCompareBin); 

						if( vpairSparse.end() ==  iterLowerBound ) 
						{
							dCumsum = vpairSparse[vpairSparse.size() - 1].second;
						} 
						else
						if( vpairSparse.begin() != iterLowerBound )
						{
							size_t uOffset = (size_t)(iterLowerBound - vpairSparse.begin());
							dCumsum = vpairSparse[uOffset - 1].second;
						}

						udValues.push_back(dCumsum);
					}
				}
			}	
		}
		#endif			//	#if	!WITHOUT_BIN_AGGREGATION

		virtual	
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
				#if	WITHOUT_FULL_ARRAYS		
				CDecodingSparseArrays::iterator iterSparseArrays = pcDecodingSparseArrays->find(uIndex);
				if( pcDecodingSparseArrays->end() 
						 != iterSparseArrays && 
					NULL != iterSparseArrays->second )
				{
					const vector< pair<BT, WT> >& vpairSparse = *iterSparseArrays->second;
					vpairCoefs.resize(vpairSparse.size());
					copy(vpairSparse.begin(), vpairSparse.end(), vpairCoefs.begin());
				}
				#else	// #if	WITHOUT_FULL_ARRAYS	
				if( vcDecodingSparseArrays[uIndex] ) 
				{
					const vector< pair<BT, WT> >& vpairSparse = *vcDecodingSparseArrays[uIndex];
					vpairCoefs.resize(vpairSparse.size());
					copy(vpairSparse.begin(), vpairSparse.end(), vpairCoefs.begin());
				}
				#endif	// #if	WITHOUT_FULL_ARRAYS	
			}	

			#if	!WITHOUT_BIN_AGGREGATION
			WT dCumsum = 0.0;
			for(vector< pair<BT, WT> >::iterator 
					iterCoefs = vpairCoefs.begin();
				iterCoefs != vpairCoefs.end();
				iterCoefs ++)
			{
				WT dTemp = iterCoefs->second;
				iterCoefs->second -= dCumsum;
				dCumsum = dTemp;
			}
			#endif	// #if	!WITHOUT_BIN_AGGREGATION
		}
		
		virtual
		void
		_GetCoefSparse
		(
			const vector<size_t>& vuSubs,
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

		virtual
		~CSepDWTDecoderPool()
		{
			#if	WITHOUT_FULL_ARRAYS		
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
			#else	// #if	WITHOUT_FULL_ARRAYS	
			for(CDecodingSparseArrays::iterator 
					iterSparseArrays = vcDecodingSparseArrays.begin(); 
				iterSparseArrays != vcDecodingSparseArrays.end(); 
				iterSparseArrays++) 
			{
				if( *iterSparseArrays ) 
				{
					delete *iterSparseArrays;
				}
			}
			#endif	// #if	WITHOUT_FULL_ARRAYS	
		}
	};
}
		

/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
