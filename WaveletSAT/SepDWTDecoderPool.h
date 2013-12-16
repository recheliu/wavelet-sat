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
		#if	WITHOUT_FULL_ARRAYS	// ADD-BY-LEETEN 12/07/2013
		typedef unordered_map< size_t, CDecodingSparseArray* > CDecodingSparseArrays;
		CDecodingSparseArrays *pcDecodingSparseArrays;
		// ADD-BY-LEETEN 12/07/2013-BEGIN
		#else	// #if	WITHOUT_FULL_ARRAYS	
		typedef vector< CDecodingSparseArray* > CDecodingSparseArrays;
		CDecodingSparseArrays vcDecodingSparseArrays;
		#endif	// #if	WITHOUT_FULL_ARRAYS	
		// ADD-BY-LEETEN 12/07/2013-END

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
			#if	WITHOUT_FULL_ARRAYS		// ADD-BY-LEETEN 12/07/2013
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
			// ADD-BY-LEETEN 12/07/2013-BEGIN
			#else	// #if	WITHOUT_FULL_ARRAYS	
			CDecodingSparseArray* pcSparseArray = vcDecodingSparseArrays[uIndex];
			if( !pcSparseArray )
			{
				pcSparseArray = new CDecodingSparseArray();
				vcDecodingSparseArrays[uIndex] = pcSparseArray;
			}
			#endif	// #if	WITHOUT_FULL_ARRAYS	
			// ADD-BY-LEETEN 12/07/2013-END
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
			#if	WITHOUT_FULL_ARRAYS		// ADD-BY-LEETEN 12/07/2013
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
			// ADD-BY-LEETEN 12/07/2013-BEGIN
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
			// ADD-BY-LEETEN 12/07/2013-END
		}

	public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/
		// ADD-BY-LEETEN 2013/11/02/2013-BEGIN
		void
		_Finalize
		(
			WT WaveletWeight,
			void* _Reserved = NULL
		)
		{
			if( bIsSparse ) {
				#if	WITHOUT_FULL_ARRAYS		// ADD-BY-LEETEN 12/07/2013
				for(CDecodingSparseArrays::iterator 
						iterSparseArrays = pcDecodingSparseArrays->begin(); 
					iterSparseArrays != pcDecodingSparseArrays->end(); 
					iterSparseArrays++) 
				{
					if( iterSparseArrays->second ) 
					{
						CDecodingSparseArray& sparseArray = *iterSparseArrays->second;
				// ADD-BY-LEETEN 12/07/2013-BEGIN
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
				// ADD-BY-LEETEN 12/07/2013-END
						sort(sparseArray.begin(), sparseArray.end());
		
						// ADD-BY-LEETEN 2013/12/01-BEGIN
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
						// ADD-BY-LEETEN 2013/12/01-END
					}
				}
			}

			// ADD-BY-LEETEN 2013/12/01-BEGIN
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
			// ADD-BY-LEETEN 2013/12/01-END
		}
		// ADD-BY-LEETEN 2013/11/02/2013-END

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
				#if	WITHOUT_FULL_ARRAYS		// ADD-BY-LEETEN 12/07/2013
				this->pcDecodingSparseArrays = new CDecodingSparseArrays();
				// ADD-BY-LEETEN 12/07/2013-BEGIN
				#else	// #if	WITHOUT_FULL_ARRAYS	
				vcDecodingSparseArrays.resize(UGetProduct(vuLengths));
				#endif	// #if	WITHOUT_FULL_ARRAYS	
				// ADD-BY-LEETEN 12/07/2013-END
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

		// ADD-BY-LEETEN 2013/12/01-BEGIN
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
						WT dCumsum = ( !usBin ) ? 0.0 : vFullArray[min(usBin - 1, vFullArray.size()-1)];
						udValues.push_back(dCumsum);
					}
				}
			}
			else
			{
				#if	WITHOUT_FULL_ARRAYS		// ADD-BY-LEETEN 12/07/2013
				CDecodingSparseArrays::iterator iterSparseArrays = pcDecodingSparseArrays->find(uIndex);
				if( pcDecodingSparseArrays->end() 
						 != iterSparseArrays && 
					NULL != iterSparseArrays->second )
				{
					const vector< pair<BT, WT> >& vpairSparse = *iterSparseArrays->second;
				// ADD-BY-LEETEN 12/07/2013-BEGIN
				#else	// #if	WITHOUT_FULL_ARRAYS	
				if( vcDecodingSparseArrays[uIndex] )
				{
					const vector< pair<BT, WT> >& vpairSparse = *vcDecodingSparseArrays[uIndex];
				#endif	// #if	WITHOUT_FULL_ARRAYS	
				// ADD-BY-LEETEN 12/07/2013-END

					for(vector<BT>::const_iterator 
							iterBins = vsBins.begin();
						iterBins != vsBins.end();
						iterBins ++ ) 
					{
						BT usBin = *iterBins;
						WT dCumsum = 0.0;
						#if	0	// MOD-BY-LEETEN 2013/12/14-FROM:
						if( usBin ) 
						{
							usBin--;
							// Among the bins smaller than usBin, find the largest one.
							pair<BT, WT> pairValue = make_pair<BT, WT>(usBin, (WT)0.0);	
							vector< pair<BT, WT> >::const_iterator iterLowerBound = 
								lower_bound(vpairSparse.begin(), vpairSparse.end(), pairValue, BCompareBin); 

							dCumsum = ( iterLowerBound != vpairSparse.end() )?
								iterLowerBound->second:
								vpairSparse[vpairSparse.size() - 1].second;
						}
						#else	// MOD-BY-LEETEN 2013/12/14-TO:
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
						#endif	// MOD-BY-LEETEN 2013/12/14-END

						udValues.push_back(dCumsum);
					}
				}
			}	
		}
		#endif			//	#if	!WITHOUT_BIN_AGGREGATION
		// ADD-BY-LEETEN 2013/12/01-END

		virtual	// ADD-BY-LEETEN 2013/10/30
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
				#if	WITHOUT_FULL_ARRAYS		// ADD-BY-LEETEN 12/07/2013
				CDecodingSparseArrays::iterator iterSparseArrays = pcDecodingSparseArrays->find(uIndex);
				if( pcDecodingSparseArrays->end() 
						 != iterSparseArrays && 
					NULL != iterSparseArrays->second )
				{
					const vector< pair<BT, WT> >& vpairSparse = *iterSparseArrays->second;
					vpairCoefs.resize(vpairSparse.size());
					copy(vpairSparse.begin(), vpairSparse.end(), vpairCoefs.begin());
				}
				// ADD-BY-LEETEN 12/07/2013-BEGIN
				#else	// #if	WITHOUT_FULL_ARRAYS	
				if( vcDecodingSparseArrays[uIndex] ) 
				{
					const vector< pair<BT, WT> >& vpairSparse = *vcDecodingSparseArrays[uIndex];
					vpairCoefs.resize(vpairSparse.size());
					copy(vpairSparse.begin(), vpairSparse.end(), vpairCoefs.begin());
				}
				#endif	// #if	WITHOUT_FULL_ARRAYS	
				// ADD-BY-LEETEN 12/07/2013-END
			}	

			// ADD-BY-LEETEN 2013/12/01-BEGIN
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
			// ADD-BY-LEETEN 2013/12/01-END
		}
		
		virtual
		void
		_GetCoefSparse
		(
			// MOD-BY-LEETEN 2013/10/30:			const vector<size_t> vuSubs,
			const vector<size_t>& vuSubs,
			// MOD-BY-LEETEN 2013/10/30-END
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
			#if	WITHOUT_FULL_ARRAYS		// ADD-BY-LEETEN 12/07/2013
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
			// ADD-BY-LEETEN 12/07/2013-BEGIN
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
			// ADD-BY-LEETEN 12/07/2013-END
		}
		// ADD-BY-LEETEN 11/12/2012-END
	};
}
		
