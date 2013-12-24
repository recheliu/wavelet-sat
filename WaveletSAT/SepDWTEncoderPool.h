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
	class CSepDWTEncoderPool
		:virtual public CSepDWTPoolBase<WT, BT>
	{
	protected:	
		// In the pair, the first element record the current count
		typedef pair<size_t, unordered_map<BT, WT>*> CEncodingSparseArray;
		typedef unordered_map<size_t, CEncodingSparseArray*> CEncodingSparseArrays; 
		CEncodingSparseArrays *pcEncodingSparseArrays;
		
		// ADD-BY-LEETEN 2013/07/12-BEGIN
		#if	WITH_STREAMING
		struct CFileBuffer {
			vector<size_t> vuHeaderGlobalOffsets;
			vector<size_t> vuHeaderLengthsNotInBuffer;
			vector<size_t> vuHeaderLengthsInBuffer;

			size_t uNrOfBufferedHeaders;
			vector<CSATSepDWTNetCDF::TYPE_COEF_COUNT>		vCounts;
			vector<CSATSepDWTNetCDF::TYPE_COEF_OFFSET>	vOffsets;
			vector<CSATSepDWTNetCDF::TYPE_COEF_BIN>		vBins;
			vector<CSATSepDWTNetCDF::TYPE_COEF_VALUE>		vValues;

			size_t uNrOfFlushes;
			size_t uMaxNrOfHeaders;
		} cFileBuffer;
		#endif	// #if	WITH_STREAMING
		// ADD-BY-LEETEN 2013/07/12-END

	protected:
		void
		_MoveToBuffer(
			const size_t& uIndex,
			void* _Reserved = NULL
		)
		{
			#if	WITH_STREAMING
			size_t uCount = 0;
			cFileBuffer.vOffsets[cFileBuffer.uNrOfBufferedHeaders] = cFileBuffer.vValues.size();
			#endif	//	#if	WITH_STREAMING
			if( !bIsSparse ) 
			{
				CFullArrays::iterator iterFullArrays = this->pcFullArrays->find(uIndex);
				if( pcFullArrays->end()
						 != iterFullArrays &&
					NULL != iterFullArrays->second &&
					NULL != iterFullArrays->second->second )
				{
					vector<WT>& vFullArray = *iterFullArrays->second->second;
					for(size_t b = 0; b < vFullArray.size(); b++)
					{
						if( vFullArray[b] )
						{
						vFullArray[b] *= this->dWaveletWeight;
					#if	WITH_STREAMING
							cFileBuffer.vBins.push_back(	(CSATSepDWTNetCDF::TYPE_COEF_BIN)b );
							cFileBuffer.vValues.push_back(	(CSATSepDWTNetCDF::TYPE_COEF_VALUE)vFullArray[b] );
							uCount++;
					#endif	// #if	WITH_STREAMING
						}
					}
					
					#if	WITH_STREAMING
					delete iterFullArrays->second->second;
					delete iterFullArrays->second;
					pcFullArrays->erase(iterFullArrays);
					#endif	// #if	WITH_STREAMING
				}
			}
			else
			{
				CEncodingSparseArrays::iterator iterSparseArrays = this->pcEncodingSparseArrays->find(uIndex);
				if( pcEncodingSparseArrays->end()
						 != iterSparseArrays &&
					NULL != iterSparseArrays->second &&
					NULL != iterSparseArrays->second->second )
				{
					unordered_map<BT, WT>& mapSparseArray = *iterSparseArrays->second->second;
					for(unordered_map<BT, WT>::iterator
							iterSparseArray = mapSparseArray.begin();
						iterSparseArray != mapSparseArray.end();
						iterSparseArray++)
					{
						iterSparseArray->second *= dWaveletWeight;
						#if	WITH_STREAMING
						cFileBuffer.vBins.push_back(	(CSATSepDWTNetCDF::TYPE_COEF_BIN)iterSparseArray->first );
						cFileBuffer.vValues.push_back(	(CSATSepDWTNetCDF::TYPE_COEF_VALUE)iterSparseArray->second );
						uCount++;
						#endif	// #if	WITH_STREAMING
					}
					#if	WITH_STREAMING
					delete iterSparseArrays->second->second;
					delete iterSparseArrays->second;
					pcEncodingSparseArrays->erase(iterSparseArrays);
					#endif	// #if	WITH_STREAMING
				}
			}
			#if	WITH_STREAMING
			cFileBuffer.vCounts[cFileBuffer.uNrOfBufferedHeaders] = (CSATSepDWTNetCDF::TYPE_COEF_COUNT)uCount;
			cFileBuffer.uNrOfBufferedHeaders++;
			#endif	// #if	WITH_STREAMING
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
			CSepDWTPoolBase::_AddAtFullArray(
				usBin,
				uIndex,
				Value,
				uCount,
				_Reserved
			);

			CFullArrays::iterator iterFullArrays = this->pcFullArrays->find(uIndex);
			CFullArray* pcFullArray = pcFullArray = iterFullArrays->second;
			size_t uCurrentCount = pcFullArray->first;
			if( uCurrentCount >= uMaxCount ) 
			{
				this->_MoveToBuffer(uIndex);
			}
		}

		void
		_AddAtEncodingSparseArray
		(
			const BT& usBin,
			const size_t& uIndex,
			const WT& Value,
			const size_t& uCount,	
			void* _Reserved = NULL
		) 
		{
			CEncodingSparseArrays::iterator iterSparseArrays = this->pcEncodingSparseArrays->find(uIndex);
			CEncodingSparseArray* pcSparseArray = NULL;
			if( iterSparseArrays != pcEncodingSparseArrays->end() )
			{
				pcSparseArray = iterSparseArrays->second;
			}
			else
			{
				pcSparseArray = new CEncodingSparseArray();
				this->pcEncodingSparseArrays->insert(pair<size_t, CEncodingSparseArray*>(uIndex, pcSparseArray));
			}

			if( !pcSparseArray->second ) 
			{
				pcSparseArray->first = UGetNonOccupiedVolume(uIndex);
				pcSparseArray->second = new unordered_map<BT, WT>();
			} 
			pcSparseArray->first += uCount;
			unordered_map<BT, WT>& mapSparseArray = *pcSparseArray->second;
			unordered_map<BT, WT>::iterator iterSparseArray = mapSparseArray.find(usBin);
			if( iterSparseArray != mapSparseArray.end() ) 
			{
				iterSparseArray->second += Value; 
			}
			else
			{
				if( Value )	// ADD-BY-LEETEN 2013/07/08
					mapSparseArray.insert(pair<BT, WT>(usBin, Value));
			}

			if( !uCount ) 
				return;

			size_t uCurrentCount = pcSparseArray->first;
			if( uCurrentCount >= uMaxCount ) 
			{
				this->_MoveToBuffer(uIndex);
			}
		}

		void
		_GetAtEncodingSparseArray
		(
			const BT& Bin,
			const size_t& uIndex,
			WT& Value,
			void* _Reserved = NULL
		) 
		{
			CEncodingSparseArrays::iterator iterEncodingSparseArrays = pcEncodingSparseArrays->find(uIndex);
			if( 
				pcEncodingSparseArrays->end()
					 == iterEncodingSparseArrays ||
				NULL == iterEncodingSparseArrays->second ||
				NULL == iterEncodingSparseArrays->second->second )
			{
				Value = WT(0);
				return;
			}

			unordered_map<BT, WT>& mapSparseArray = *iterEncodingSparseArrays->second->second;
			unordered_map<BT, WT>::iterator iterSparseArray = mapSparseArray.find(Bin);
			Value = ( mapSparseArray.end() == iterSparseArray )?WT(0):iterSparseArray->second;
		}

	public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		#if	WITH_STREAMING
		void
		_SetBuffer
		(
			size_t uMinNrOFBufferedThreads,
			const vector<size_t>& vuHeaderGlobalOffsets,
			void* _Reserved = NULL
		)
		{
			this->cFileBuffer.vuHeaderGlobalOffsets.assign(vuHeaderGlobalOffsets.begin(), vuHeaderGlobalOffsets.end());

			cFileBuffer.uNrOfBufferedHeaders = 0;
			cFileBuffer.uNrOfFlushes = 0;

			cFileBuffer.vuHeaderLengthsNotInBuffer.resize(vuDataDimLengths.size());
			cFileBuffer.vuHeaderLengthsInBuffer.resize(vuDataDimLengths.size());
			cFileBuffer.uMaxNrOfHeaders = 1;
			for(size_t d = 0; d < vuDataDimLengths.size(); d++) 
			{
				size_t uHeaderLength = (size_t)ceilf((float)vuDataDimLengths[d] / (float)vuWaveletLengths[d]);
				if( uMinNrOFBufferedThreads >= cFileBuffer.uMaxNrOfHeaders )
				{
					cFileBuffer.uMaxNrOfHeaders *= uHeaderLength;
					cFileBuffer.vuHeaderLengthsInBuffer[d] = uHeaderLength;
					cFileBuffer.vuHeaderLengthsNotInBuffer[d] = 1;
				} 
				else 
				{
					cFileBuffer.vuHeaderLengthsInBuffer[d] = 1;
					cFileBuffer.vuHeaderLengthsNotInBuffer[d] = uHeaderLength;
				}
			}

			cFileBuffer.vCounts.resize(cFileBuffer.uMaxNrOfHeaders);
			cFileBuffer.vOffsets.resize(cFileBuffer.uMaxNrOfHeaders);
		}

		bool 
		BIsReadyToFlush
		(
			void *_Reserved = NULL
		)
		{
			if( cFileBuffer.uNrOfBufferedHeaders == cFileBuffer.uMaxNrOfHeaders )
				return true;
			return false;
		}

		void
		_GetFileBuffer(
			size_t puHeaderStart[NC_MAX_DIMS],
			size_t puHeaderCount[NC_MAX_DIMS], 
			size_t& uNrOfBufferedHeaders,
			CSATSepDWTNetCDF::TYPE_COEF_COUNT	**ppCoefCounts,
			CSATSepDWTNetCDF::TYPE_COEF_OFFSET	**ppCoefOffsets,

			size_t& uNrOfBufferedCoefs,
			CSATSepDWTNetCDF::TYPE_COEF_VALUE	**ppCoefValues,
			CSATSepDWTNetCDF::TYPE_COEF_BIN		**ppCoefBins,

			void *_Reserved = NULL
			)
		{
			// convert the current scanline to its indices
			vector<size_t> vuFlushOffsets;
			_ConvertIndexToSub(cFileBuffer.uNrOfFlushes, vuFlushOffsets, cFileBuffer.vuHeaderLengthsNotInBuffer);
			size_t uNrOfDims = cFileBuffer.vuHeaderGlobalOffsets.size();
			for(size_t d = 0; d < uNrOfDims; d++) {
			    puHeaderStart[uNrOfDims - 1 - d] = cFileBuffer.vuHeaderGlobalOffsets[d] + vuFlushOffsets[d];
			    puHeaderCount[uNrOfDims - 1 - d] = cFileBuffer.vuHeaderLengthsInBuffer[d];
			}

			uNrOfBufferedHeaders = cFileBuffer.uNrOfBufferedHeaders;
			uNrOfBufferedCoefs = cFileBuffer.vValues.size();
			*ppCoefCounts	= cFileBuffer.vCounts.data();
			*ppCoefOffsets	= cFileBuffer.vOffsets.data();
			*ppCoefValues	= cFileBuffer.vValues.data();
			*ppCoefBins		= cFileBuffer.vBins.data();
			cFileBuffer.uNrOfFlushes++;
		}

		void
		_ResetFileBuffer(
			void *_Reserved = NULL
		)
		{
			cFileBuffer.uNrOfBufferedHeaders = 0;
			cFileBuffer.vValues.clear();
			cFileBuffer.vBins.clear();
		}
		#endif	// #if	WITH_STREAMING

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
				bIsSparse
			);

			if( bIsSparse )
			{	
				this->pcEncodingSparseArrays = new CEncodingSparseArrays();
			}
		}

		void
		_GetArraySize
		(
			size_t& uCountInFullArray,
			size_t& uCountInSparseArray,
			const WT& Threshold,
			void* _Reserved = NULL
		)
		{
			uCountInFullArray = 0;
			uCountInSparseArray = 0;
			for(size_t i = 0; i < uSize; i++) 
			{
				if(!bIsSparse)
				{
					vector<pair<BT, WT>> vpair;
					this->_GetCoefSparse(i, vpair);
					uCountInFullArray  += vpair.size();
				}
				else
				{
					CEncodingSparseArrays::iterator iterSparseArrays = pcEncodingSparseArrays->find(i);
					if( pcEncodingSparseArrays->end() 
							 != iterSparseArrays && 
						NULL != iterSparseArrays->second &&
						NULL != iterSparseArrays->second->second ) 
					{
						uCountInSparseArray  += iterSparseArrays->second->second->size();
					}
				}
			}
		}

		void
		_GetAt
		(
			const BT& Bin,
			const vector<size_t>& vuSubs,
			size_t& uIndex,
			WT& Value,
			void* _Reserved = NULL
		)
		{
			if( vuSubs.size() )
				uIndex = UConvertSubToIndex(vuSubs, vuLengths);

			if( !bIsSparse )
			{
				_GetAtFullArray(Bin, uIndex, Value);
			}
			else
			{
				_GetAtEncodingSparseArray(Bin, uIndex, Value);
			}
		}

		//! Add value to the location specified by the 1D index
		void
		_AddAt
		(
			const BT& Bin,
			const vector<size_t>& vuSubs,
			const WT& Value,
			const size_t& uCount = 1,	
			void* _Reserved = NULL
		)
		{
			size_t uIndex = UConvertSubToIndex(vuSubs, vuLengths);

			if( !bIsSparse )
				this->_AddAtFullArray(Bin, uIndex, Value, uCount);
			else
				this->_AddAtEncodingSparseArray(Bin, uIndex, Value, uCount);
		}

		const 
		bool
		BIsEmpty(
			size_t uIndex,
			void* _Reserved = NULL
		) const
		{
			if( !bIsSparse )
			{
				if(!pcFullArrays)
					return true;
				CFullArrays::iterator iterFullArrays = pcFullArrays->find(uIndex);
				if(pcFullArrays->end() 
						 != iterFullArrays &&
					NULL != iterFullArrays->second )
					return (!iterFullArrays->second->first)?true:false;
				return false;
			}
			else
			{
				if(!this->pcEncodingSparseArrays )
					return true;
				CEncodingSparseArrays::iterator iterSparseArrays = pcEncodingSparseArrays->find(uIndex);
				if(pcEncodingSparseArrays->end() 
						 != iterSparseArrays &&
					NULL != iterSparseArrays->second )
					return (!iterSparseArrays->second->first)?true:false;
				return false;
			}
		}

		const 
		bool
		BIsEmpty(
			const vector<size_t> vuSubs,
			void* _Reserved = NULL
		) const
		{
			return BIsEmpty(UConvertSubToIndex(vuSubs, vuLengths));
		}

		CSepDWTEncoderPool():
			CSepDWTPoolBase()
		{
			this->pcEncodingSparseArrays = NULL;
		}

		virtual	
		~CSepDWTEncoderPool()
		{
			if( this->pcEncodingSparseArrays ) {
				for(CEncodingSparseArrays::iterator 
						iterSparseArrays = pcEncodingSparseArrays->begin(); 
					iterSparseArrays != pcEncodingSparseArrays->end(); 
					iterSparseArrays++) 
				{
					if( iterSparseArrays->second ) 
					{
						if( iterSparseArrays->second->second ) 
							delete iterSparseArrays->second->second;
						delete iterSparseArrays->second;
					}
				}
			}
		}
	};
}
		

/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
