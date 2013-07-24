#pragma once

#include	<algorithm>	// ADD-BY-LEETEN 03/17/2013

#include <unordered_map>	
#include <vector>
using namespace std;
#include <math.h>

#include "HeaderBase.h"	
#include "DWT.h"
#include "SATSepDWTNetCDF.h"	// ADD-BY-LEETEN 01/27/2013
#include "liblog.h"	

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

	//! SepDWT coefficients per basis of all bins
	template<
		typename WT = typeWavelet,	//!< Type of the wavelet coefficientsd
		typename BT = typeBin		//!< Type of the bin
	>
	class CSepDWTPool
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

		// ADD-BY-LEETEN 2013/07/08-BEGIN
		vector<size_t> vuDataDimLengths;

		vector<size_t> vuWaveletLengths;
		// ADD-BY-LEETEN 2013/07/08-END

		#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:	#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
		vector< vector<WT> > vvFull;

		vector< unordered_map<BT, WT>* >* pvpmapSparse;
		#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:			#else	// MOD-BY-LEETEN 2013/07/07-TO:
		WT	 dWaveletWeight;

		size_t uNrOfBins;

		typedef pair<size_t, vector<WT>*> CFullArray;
		typedef unordered_map<size_t,  CFullArray*> CFullArrays;
		CFullArrays* pcFullArrays;

		// In the pair, the first element record the current count
		typedef pair<size_t, unordered_map<BT, WT>*> CEncodingSparseArray;
		typedef unordered_map<size_t, CEncodingSparseArray*> CEncodingSparseArrays; 
		CEncodingSparseArrays *pcEncodingSparseArrays;
		#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:			#endif	// MOD-BY-LEETEN 2013/07/07-END
		
		// ADD-BY-LEETEN 11/11/2012-BEGIN
		//! Max # that each coefficient is updated
		size_t uMaxCount;

		#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:		#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
		//! Record how many time each coefficient has been updated
		vector<size_t> vuCounts;

		//! 
		vector< vector< pair<BT, WT> > > vvpairSparse;
		// ADD-BY-LEETEN 11/11/2012-END
		#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:		#else	// MOD-BY-LEETEN 2013/07/07-TO:
		typedef vector<pair<BT, WT>> CDecodingSparseArray;
		typedef unordered_map< size_t, CDecodingSparseArray* > CDecodingSparseArrays;
		CDecodingSparseArrays *pcDecodingSparseArrays;

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

		void
		// MOD-BY-LEETEN 2013/07/12-FROM:	_MoveToPool(
		_MoveToBuffer(
		// MOD-BY-LEETEN 2013/07/12-END
			const size_t& uIndex,
			void* _Reserved = NULL
		)
		{
			// ADD-BY-LEETEN 2013/07/12-BEGIN
			#if	WITH_STREAMING
			size_t uCount = 0;
			cFileBuffer.vOffsets[cFileBuffer.uNrOfBufferedHeaders] = cFileBuffer.vValues.size();
			#endif	//	#if	WITH_STREAMING
			// ADD-BY-LEETEN 2013/07/12-END
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
					// ADD-BY-LEETEN 2013/07/12-BEGIN
					{
						if( vFullArray[b] )
						{
					// ADD-BY-LEETEN 2013/07/12-END
						vFullArray[b] *= this->dWaveletWeight;
					// ADD-BY-LEETEN 2013/07/12-BEGIN
					#if	WITH_STREAMING
							cFileBuffer.vBins.push_back(	(CSATSepDWTNetCDF::TYPE_COEF_BIN)b );
							cFileBuffer.vValues.push_back(	(CSATSepDWTNetCDF::TYPE_COEF_VALUE)vFullArray[b] );
							uCount++;
					#endif	// #if	WITH_STREAMING
						}
					}
					// ADD-BY-LEETEN 2013/07/12-END
					
					#if	0	// MOD-BY-LEETEN 2013/07/12-FROM:
					// TODO: Move to the pool
					/*
					delete iterFullArrays->second->second;
					delete iterFullArrays->second;
					pcFullArrays->erase(iterFullArrays);
					*/
					#else	// MOD-BY-LEETEN 2013/07/12-TO:
					#if	WITH_STREAMING
					delete iterFullArrays->second->second;
					delete iterFullArrays->second;
					pcFullArrays->erase(iterFullArrays);
					#endif	// #if	WITH_STREAMING
					#endif	// MOD-BY-LEETEN 2013/07/12-END
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
						// ADD-BY-LEETEN 2013/07/12-BEGIN
						#if	WITH_STREAMING
						cFileBuffer.vBins.push_back(	(CSATSepDWTNetCDF::TYPE_COEF_BIN)iterSparseArray->first );
						cFileBuffer.vValues.push_back(	(CSATSepDWTNetCDF::TYPE_COEF_VALUE)iterSparseArray->second );
						uCount++;
						#endif	// #if	WITH_STREAMING
						// ADD-BY-LEETEN 2013/07/12-END
					}
					#if	0	// MOD-BY-LEETEN 2013/07/12-FROM:
					// TODO: Move to the pool
					/*
					delete iterSparseArrays->second->second;
					delete iterSparseArrays->second;
					pcEncodingSparseArrays->erase(iterSparseArrays);
					*/
					#else	// MOD-BY-LEETEN 2013/07/12-TO:
					#if	WITH_STREAMING
					delete iterSparseArrays->second->second;
					delete iterSparseArrays->second;
					pcEncodingSparseArrays->erase(iterSparseArrays);
					#endif	// #if	WITH_STREAMING
					#endif	// MOD-BY-LEETEN 2013/07/12-END
				}
			}
			// ADD-BY-LEETEN 2013/07/12-BEGIN
			#if	WITH_STREAMING
			cFileBuffer.vCounts[cFileBuffer.uNrOfBufferedHeaders] = (CSATSepDWTNetCDF::TYPE_COEF_COUNT)uCount;
			cFileBuffer.uNrOfBufferedHeaders++;
			#endif	// #if	WITH_STREAMING
			// ADD-BY-LEETEN 2013/07/12-END
		}

		// ADD-BY-LEETEN 2013/07/08-BEGIN
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
		// ADD-BY-LEETEN 2013/07/08-END

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
				// MOD-BY-LEETEN 2013/07/08-FROM:				pcFullArray->first = 0;
				pcFullArray->first = UGetNonOccupiedVolume(uIndex);
				// MOD-BY-LEETEN 2013/07/08-END

				pcFullArray->second = new vector<WT>();
				pcFullArray->second->assign(uNrOfBins, (WT)0);
			} 
			pcFullArray->first += uCount;
			vector<WT>& vec = *pcFullArray->second;
			vec[usBin] += Value;

			if( !uCount ) 
				return;

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
				// MOD-BY-LEETEN 2013/07/08-FROM:				pcSparseArray->first = 0;
				pcSparseArray->first = UGetNonOccupiedVolume(uIndex);
				// MOD-BY-LEETEN 2013/07/08-END
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
		#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:		#endif	// MOD-BY-LEETEN 2013/07/07-END

	public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0800,
			PARAMETER_END
		};

		bool 
		BIsSparse
		(
			void* _Reserved = NULL
		)
		{
			return bIsSparse;
		}
		
		// ADD-BY-LEETEN 2013/07/12-BEGIN
		#if	WITH_STREAMING
		void
		_SetBuffer
		(
			// ADD-BY-LEETEN 2013/07/14-BEGIN
			size_t uMinNrOFBufferedThreads,
			// ADD-BY-LEETEN 2013/07/14-END
			const WT& dWaveletWeight,
			const vector<size_t>& vuDataDimLengths,
			const vector<size_t>& vuWaveletLengths,
			const vector<size_t>& vuHeaderGlobalOffsets,
			void* _Reserved = NULL
		)
		{
			this->dWaveletWeight = dWaveletWeight;	
			this->vuDataDimLengths.assign(vuDataDimLengths.begin(), vuDataDimLengths.end());
			this->vuWaveletLengths.assign(vuWaveletLengths.begin(), vuWaveletLengths.end());
			this->cFileBuffer.vuHeaderGlobalOffsets.assign(vuHeaderGlobalOffsets.begin(), vuHeaderGlobalOffsets.end());

			cFileBuffer.uNrOfBufferedHeaders = 0;
			cFileBuffer.uNrOfFlushes = 0;

			cFileBuffer.vuHeaderLengthsNotInBuffer.resize(vuDataDimLengths.size());
			cFileBuffer.vuHeaderLengthsInBuffer.resize(vuDataDimLengths.size());
			cFileBuffer.uMaxNrOfHeaders = 1;
			for(size_t d = 0; d < vuDataDimLengths.size(); d++) 
			{
				size_t uHeaderLength = (size_t)ceilf((float)vuDataDimLengths[d] / (float)vuWaveletLengths[d]);
				#if	0	// MOD-BY-LEETEN 2013/07/14-FROM:
				#if	1	// TMP-MOD
				if( 1 == cFileBuffer.uMaxNrOfHeaders )
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
				#else
				cFileBuffer.uMaxNrOfHeaders *= uHeaderLength;
				cFileBuffer.vuHeaderLengthsInBuffer[d] = uHeaderLength;
				cFileBuffer.vuHeaderLengthsNotInBuffer[d] = 1;
				#endif
				#else	// MOD-BY-LEETEN 2013/07/14-TO:
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
				#endif	// MOD-BY-LEETEN 2013/07/14-END
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
		// ADD-BY-LEETEN 2013/07/12-END

		// ADD-BY-LEETEN 01/27/2013-BEGIN
		void
		_Copy(
			size_t uIndex,
			size_t uNrOfBins,
			const CSATSepDWTNetCDF::TYPE_COEF_BIN	*pBins,
			const CSATSepDWTNetCDF::TYPE_COEF_VALUE	*pValues,
			void *_Reserved = NULL
		)
		{
			#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
			if( bIsSparse )
				vvpairSparse[uIndex].resize(uNrOfBins);

			for(size_t b = 0; b < uNrOfBins; b++)
			{
				CSATSepDWTNetCDF::TYPE_COEF_VALUE	Value = pValues[b];
				CSATSepDWTNetCDF::TYPE_COEF_BIN		Bin =	pBins[b];
				if( Value )
				{
					if( !bIsSparse )
						vvFull[Bin][uIndex] = Value;
					else
						vvpairSparse[uIndex][b] = pair<BT, WT>((BT)Bin, (WT)Value);
				}	
			}

			if( bIsSparse )
				vuCounts[uIndex] = uMaxCount;
			#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#else	// MOD-BY-LEETEN 2013/07/07-TO:
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
			#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#endif	// MOD-BY-LEETEN 2013/07/07-END
		}
		// ADD-BY-LEETEN 01/27/2013-END

		#if	WITH_DYNAMIC_ARRAY_ALLOCATION		// ADD-BY-LEETEN 2013/07/23
		// ADD-BY-LEETEN 2013/07/08-BEGIN
		#if	!WITH_STREAMING	// ADD-BY-LEETEN 2013/07/12
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
		#endif	// #if	!WITH_STREAMING	// ADD-BY-LEETEN 2013/07/12
		// ADD-BY-LEETEN 2013/07/08-END
		#endif	// #if	WITH_DYNAMIC_ARRAY_ALLOCATION		// ADD-BY-LEETEN 2013/07/23

		void
		_Set
		(
			const BT& uNrOfBins,
			const vector<size_t>& vuLengths,
			size_t uMaxCount,	// ADD-BY-LEETEN 11/11/2012
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
			#if	WITH_DYNAMIC_ARRAY_ALLOCATION		// ADD-BY-LEETEN 2013/07/23
			this->uNrOfBins = uNrOfBins;	// ADD-BY-LEETEN 2013/07/07
			#endif	// #if	WITH_DYNAMIC_ARRAY_ALLOCATION		// ADD-BY-LEETEN 2013/07/23
			if( !bIsSparse )
			{
				#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:					#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
				vvFull.resize(uNrOfBins);
				for(size_t b = 0; b < (size_t)uNrOfBins; b++)
					vvFull[b].resize(this->uSize);
				#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:					#else	// MOD-BY-LEETEN 2013/07/07-TO:
				this->pcFullArrays = new CFullArrays();
				#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:					#endif	// MOD-BY-LEETEN 2013/07/07-END
			}
			else
			{	// ADD-BY-LEETEN 11/11/2012
				#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:					#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
				this->pvpmapSparse = new vector< unordered_map <BT, WT>* >;
				this->pvpmapSparse->resize(this->uSize);

			// ADD-BY-LEETEN 11/11/2012-BEGIN
				this->vvpairSparse.resize(this->uSize);
				#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:					#else	// MOD-BY-LEETEN 2013/07/07-TO:
				this->pcEncodingSparseArrays = new CEncodingSparseArrays();
				this->pcDecodingSparseArrays = new CDecodingSparseArrays();
				#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:					#endif	// MOD-BY-LEETEN 2013/07/07-END
			}
			// ADD-BY-LEETEN 11/11/2012-END
			// DEL-BY-LEETEN 2013/07/07:	this->vuCounts.assign(this->uSize, 0);
			// ADD-BY-LEETEN 2013/07/23-BEGIN
			#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
			this->vuCounts.assign(this->uSize, 0);
			#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
			// ADD-BY-LEETEN 2013/07/23-END

			this->uMaxCount = uMaxCount;	// ADD-BY-LEETEN 11/11/2012
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
			#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
			if(!bIsSparse)
			{
				for(typename vector< vector<WT> >::iterator 
					ivvFull = this->vvFull.begin();
					ivvFull != this->vvFull.end();
					ivvFull++)
				{
					vector<WT>& vBinFull = *ivvFull;
					for(typename vector<WT>::iterator 
						ivFull = vBinFull.begin();
						ivFull != vBinFull.end();
						ivFull++)
						if( (WT)fabs((double)*ivFull) >= (WT)Threshold )
							uCountInFullArray++;
				}
			}
			else
			{
				for(typename vector< vector< pair<BT, WT> > >::iterator 
					ivvpairSparse = this->vvpairSparse.begin();
					ivvpairSparse != this->vvpairSparse.end();
					ivvpairSparse++)
				{
					vector< pair<BT, WT> >& vpairSparse = *ivvpairSparse;
					for(typename vector< pair<BT, WT> >::iterator 
						pair = vpairSparse.begin();
						pair != vpairSparse.end();
						pair++)
						if( (WT)fabs((double)pair->second) >= (WT)Threshold )
							uCountInSparseArray++;
				}
				// ADD-BY-LEETEN 11/11/2012-END
			}
			#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#else	// MOD-BY-LEETEN 2013/07/07-TO:
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
			#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#endif	// MOD-BY-LEETEN 2013/07/07-END
		}


		#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:			#if	0		// DEL-BY-LEETEN 2013/07/07-BEGIN
		void
		_AddEntryToSparseArray
		(
			unordered_map<BT, WT>& mapSparse,
			const BT& Bin,
			const WT& Value,
			void* _Reserved = NULL
		)
		{
			if( Value )
			{
				mapSparse.insert(pair<BT, WT>(Bin, Value));
			}
		}
		#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:			#endif		// DEL-BY-LEETEN 2013/07/07-END

		// ADD-BY-LEETEN 12/29/2012-BEGIN
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
				#if	0	// MOD-BY-LEETEN 2013/07/23-FROM:
				// MOD-BY-LEETEN 2013/07/07-FROM:				Value = vvFull[Bin][uIndex];
				_GetAtFullArray(Bin, uIndex, Value);
				// MOD-BY-LEETEN 2013/07/07-END
				#else	// MOD-BY-LEETEN 2013/07/23-TO:
				#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
				Value = vvFull[Bin][uIndex];
				#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
				_GetAtFullArray(Bin, uIndex, Value);
				#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
				#endif	// MOD-BY-LEETEN 2013/07/23-END
			}
			else
			{
				// sparse
				#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
				vector< pair<BT, WT> >& vpair = this->vvpairSparse[uIndex];
				if( (size_t)usOffset < vpair.size() )
				{
					Bin = vpair[usOffset].first;
					Value = vpair[usOffset].second;
				}
				#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#else	// MOD-BY-LEETEN 2013/07/07-TO:
				_GetAtDecodingSparseArray(usOffset, uIndex, Bin, Value);
				#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#endif	// MOD-BY-LEETEN 2013/07/07-END
			}
		}
		// ADD-BY-LEETEN 12/29/2012-END

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
				// full
				#if	0	// MOD-BY-LEETEN 2013/07/23-FROM:
				// MOD-BY-LEETEN 2013/07/07-FROM:				Value = vvFull[Bin][uIndex];
				_GetAtFullArray(Bin, uIndex, Value);
				// MOD-BY-LEETEN 2013/07/07-END
				#else	// MOD-BY-LEETEN 2013/07/23-TO:
				#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
				Value = vvFull[Bin][uIndex];
				#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
				_GetAtFullArray(Bin, uIndex, Value);
				#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
				#endif	// MOD-BY-LEETEN 2013/07/23-END
			}
			else
			{
				// sparse
				#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
				if( NULL == (*pvpmapSparse)[uIndex] )
				{
					Value = WT(0);
					return;
				}

				unordered_map<BT, WT>& mapCoef = *(*pvpmapSparse)[uIndex];
				typename unordered_map<BT, WT>::iterator ipair = mapCoef.find(Bin);
				if( mapCoef.end() == ipair )
					Value = WT(0);
				else
					Value = ipair->second;
				#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:								#else	// MOD-BY-LEETEN 2013/07/07-TO:
				_GetAtEncodingSparseArray(Bin, uIndex, Value);
				#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:								#endif	// MOD-BY-LEETEN 2013/07/07-END
			}
		}

		// ADD-BY-LEETEN 04/26/2013-BEGIN
		void
		_IncreaseCount
		(
			size_t uExtraCount,
			void *_Reserved = NULL
		)
		{
			uMaxCount += uExtraCount;
		}
		// ADD-BY-LEETEN 04/26/2013-END

		//! Add value to the location specified by the 1D index
		void
		_AddAt
		(
			const BT& Bin,
			const vector<size_t>& vuSubs,
			const WT& Value,
			const size_t& uCount = 1,	// ADD-BY-LEETEN 01/13/2013
			void* _Reserved = NULL
		)
		{
			size_t uIndex = UConvertSubToIndex(vuSubs, vuLengths);
			// ADD-By-LEETEN 11/11/2012-BEGIN
			// DEL-BY-LEETEN 2013/07/07:				vuCounts[uIndex] += uCount;
			// ADD-BY-LEETEN 2013/07/23-BEGIN
			#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
			vuCounts[uIndex] += uCount;
			#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
			// ADD-BY-LEETEN 2013/07/23-END

			#if	0	// DEL-BY-LEETEN 2013/07/08-BEGIN
			if( Value )
			{
			#endif	// DEL-BY-LEETEN 2013/07/08-END
			// ADD-By-LEETEN 11/11/2012-END
			#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:			#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
			if( !bIsSparse )
			{
				// full
				vvFull[Bin][uIndex] += Value;
			}
			else
			{
				// sparse
				if(NULL == (*pvpmapSparse)[uIndex])
					(*pvpmapSparse)[uIndex] = new unordered_map<BT, WT>();
				unordered_map<BT, WT>& mapCoef = *(*pvpmapSparse)[uIndex];
				typename unordered_map<BT, WT>::iterator ipair = mapCoef.find(Bin);
				if( mapCoef.end() == ipair )
					_AddEntryToSparseArray
					(
						mapCoef,
						Bin,
						Value
					);
					else
						ipair->second += Value;

			}
			#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:						#else	// MOD-BY-LEETEN 2013/07/07-TO:
			if( !bIsSparse )
				this->_AddAtFullArray(Bin, uIndex, Value, uCount);
			else
				this->_AddAtEncodingSparseArray(Bin, uIndex, Value, uCount);
			#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:						#endif	// MOD-BY-LEETEN 2013/07/07-END
			#if	0	// DEL-BY-LEETEN 2013/07/08-BEGIN
			}	// ADD-By-LEETEN 11/11/2012
			#endif	// DEL-BY-LEETEN 2013/07/08-END
			#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:						#if	0	// DEL-BY-LEETEN 2013/07/07-BEGIN
			if( bIsSparse && vuCounts[uIndex] == uMaxCount )
			{
				// ADD-BY-LEETEN 04/21/2013-BEGIN
				if( !(*this->pvpmapSparse)[uIndex] )
					return;
				// ADD-BY-LEETEN 04/21/2013-END
				const unordered_map<BT, WT>& vmapBinSparse = *(*this->pvpmapSparse)[uIndex];

				vvpairSparse[uIndex].resize(vmapBinSparse.size());
				copy(vmapBinSparse.begin(), vmapBinSparse.end(), vvpairSparse[uIndex].begin());
				sort(vvpairSparse[uIndex].begin(), vvpairSparse[uIndex].end());
				/// now clear this map
				// _ShowMemoryUsage();
				(*this->pvpmapSparse)[uIndex]->clear();
				delete (*this->pvpmapSparse)[uIndex];
				(*this->pvpmapSparse)[uIndex] = NULL;
				// _ShowMemoryUsage();
			}	
			#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:						#endif	// DEL-BY-LEETEN 2013/07/07-END
		}

		void
		_Finalize
		(
			WT WaveletWeight,
			void* _Reserved = NULL
		)
		{
			#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:						#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
			if(!bIsSparse)
			{
				// ADD-BY-LEETEN 03/16/2013-BEGIN
				if( (WT)1.0 == WaveletWeight )
					return;
				// ADD-BY-LEETEN 03/16/2013-END
				for(typename vector< vector<WT> >::iterator 
					ivvFull = this->vvFull.begin();
					ivvFull != this->vvFull.end();
					ivvFull++)
				{
					vector<WT>& vBinFull = *ivvFull;
					for(typename vector<WT>::iterator 
						ivFull = vBinFull.begin();
						ivFull != vBinFull.end();
						ivFull++)
					{
						WT Coef = *ivFull;
						if(Coef)
							*ivFull *= WaveletWeight;
					}
				}
			}
			else
			{
				// check each element to make sure that no element is left in the map
				for(size_t e = 0; e < this->pvpmapSparse->size(); e++)
				{
					// ADD-By-LEETEN 11/11/2012-BEGIN
					// if the vector of pair is not empty, it means that the coefficients has been moved to here
					if(NULL == (*this->pvpmapSparse)[e])
						continue;
					// ADD-By-LEETEN 11/11/2012-END

					const unordered_map<BT, WT>& vmapBinSparse = *(*this->pvpmapSparse)[e];
					vvpairSparse[e].resize(vmapBinSparse.size());
					copy(vmapBinSparse.begin(), vmapBinSparse.end(), vvpairSparse[e].begin());
					sort(vvpairSparse[e].begin(), vvpairSparse[e].end());
					/// now clear this map
					(*this->pvpmapSparse)[e]->clear();
					delete (*this->pvpmapSparse)[e];
					(*this->pvpmapSparse)[e] = NULL;
				}

				// clear the vector to hold the sparse array and the counts
				this->pvpmapSparse->clear();
				delete this->pvpmapSparse;
				this->pvpmapSparse = NULL;
				this->vuCounts.clear();
				for(typename vector< vector< pair<BT, WT> > >::iterator 
					ivvpairSparse = this->vvpairSparse.begin();
					ivvpairSparse != this->vvpairSparse.end();
					ivvpairSparse++)
				{
					vector< pair<BT, WT> >& vpairSparse = *ivvpairSparse;
					for(typename vector< pair<BT, WT> >::iterator 
						pair = vpairSparse.begin();
						pair != vpairSparse.end();
						pair++)
						pair->second *= WaveletWeight;
				}
				// ADD-BY-LEETEN 11/11/2012-END
			}
			#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:						#else	// MOD-BY-LEETEN 2013/07/07-TO:
			// TODO:
			if( !bIsSparse ) 
			{	
				// TODO:
			}
			else
			{
				// TODO:
			}
			#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:						#endif	// MOD-BY-LEETEN 2013/07/07-END
		}

		// ADD-BY-LEETEN 03/16/2013-BEGIN
		const 
		bool
		BIsEmpty(
			size_t uIndex,
			void* _Reserved = NULL
		) const
		{
			// MOD-BY-LEETEN 2013/07/23-FROM: // MOD-BY-LEETEN 2013/07/07-FROM:			return (vuCounts[uIndex])?false:true;
			#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
			return (vuCounts[uIndex])?false:true;
			#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
			// MOD-BY-LEETEN 2013/07/23-END
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
			// MOD-BY-LEETEN 2013/07/07-END
			#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// ADD-BY-LEETEN 2013/07/23
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

		#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:	#if	0		// DEL-BY-LEETEN 2013/07/07-BEGIN
		virtual
		const 
		vector< pair<BT, WT> >&
		VGetCoefSparse
		(
			size_t uIndex,
			void* _Reserved = NULL
		) const
		{
			#if	WITH_DYNAMIC_ARRAY_ALLOCATION		// ADD-BY-LEETEN 2013/07/23
			// ADD-BY-LEETEN 2013/07/07-BEGIN
			static vector< pair<BT, WT> > vpairCoefs;	
			vpairCoefs.clear();
			// ADD-BY-LEETEN 2013/07/07-END
			#endif	// #if	WITH_DYNAMIC_ARRAY_ALLOCATION		// ADD-BY-LEETEN 2013/07/23
			if( !bIsSparse )
			{
				#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:					#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
				static vector< pair<BT, WT> > vpairCoefs;	
				vpairCoefs.clear();

				for(size_t b = 0; b < this->vvFull.size(); b++)
				{
					WT Coef = vvFull[b][uIndex];
					if( Coef )
						vpairCoefs.push_back(pair<BT, WT>((BT)b, Coef));
				}
				#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:					#else	// MOD-BY-LEETEN 2013/07/07-TO:
				CFullArray::iterator iterFullArrays = pcFullArrays->find(uIndex);
				if( iterFullArrays != pcFullArrays->end() && NULL != iterFullArrays->second ) 
				{
					CFullArray& vFullArray = *iterFullArrays->second;
					for(size_t b = 0; b < uNrOfBins ; b++)
					{
						WT Coef = vFullArray[b];
						if( Coef )
							vpairCoefs.push_back(pair<BT, WT>((BT)b, Coef));
					}
				}
				#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:					#endif	// MOD-BY-LEETEN 2013/07/07-END
				return vpairCoefs;
			}
			else
			#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
				return this->vvpairSparse[uIndex];
			#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#else	// MOD-BY-LEETEN 2013/07/07-TO:
			{
				CDecodingSparseArray::iterator iterSparseArrays = pcDecodingSparseArrays->find(uIndex);
				if( pcDecodingSparseArrays->end() 
						 != iterSparseArrays && 
					NULL != iterSparseArrays->second &&
					NULL != iterSparseArrays->second->second ) 
				{
					return *iterSparseArrays->second->second;
				}
			}
			return vpairCoefs;
			#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#endif	// MOD-BY-LEETEN 2013/07/07-END
		}
		
		virtual
		const 
		vector< pair<BT, WT> >&
		VGetCoefSparse
		(
			const vector<size_t> vuSubs,
			void* _Reserved = NULL
		) const
		{
			size_t uIndex = UConvertSubToIndex(vuSubs, vuLengths);
			return VGetCoefSparse(uIndex);
		}
		// ADD-BY-LEETEN 03/16/2013-END
		#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:			#endif	// DEL-BY-LEETEN 2013/07/07-END

		void
		_GetCoefSparse
		(
			size_t uIndex,
			vector< pair<BT, WT> >& vpairCoefs,
			void* _Reserved = NULL
		) const
		{
			// ADD-BY-LEETEN 11/12/2012-BEGIN
			vpairCoefs.clear();
			if( !bIsSparse )
			{
				#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:					#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
				for(size_t b = 0; b < this->vvFull.size(); b++)
				{
					WT Coef = vvFull[b][uIndex];
					if( Coef )
						vpairCoefs.push_back(pair<BT, WT>((BT)b, Coef));
				}
				#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:					#else	// MOD-BY-LEETEN 2013/07/07-TO:
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
				#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:					#endif	// MOD-BY-LEETEN 2013/07/07-END
			}
			else
			{
			// ADD-BY-LEETEN 11/12/2012-END
			// sparse
			#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#if	0		// MOD-BY-LEETEN 2013/07/07-FROM:
			const vector< pair<BT, WT> >& vpairSparse = this->vvpairSparse[uIndex];
			vpairCoefs.resize(vpairSparse.size());
			copy(vpairSparse.begin(), vpairSparse.end(), vpairCoefs.begin());
			#else // #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#else		// MOD-BY-LEETEN 2013/07/07-TO:
				CDecodingSparseArrays::iterator iterSparseArrays = pcDecodingSparseArrays->find(uIndex);
				if( pcDecodingSparseArrays->end() 
						 != iterSparseArrays && 
					NULL != iterSparseArrays->second )
				{
					const vector< pair<BT, WT> >& vpairSparse = *iterSparseArrays->second;
					vpairCoefs.resize(vpairSparse.size());
					copy(vpairSparse.begin(), vpairSparse.end(), vpairCoefs.begin());
				}
			#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#endif		// MOD-BY-LEETEN 2013/07/07-END
			}	// ADD-BY-LEETEN 11/12/2012
		}
		
		// ADD-BY-LEETEN 01/05/2013-BEGIN
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
		// ADD-BY-LEETEN 01/05/2013-END

		CSepDWTPool()
		{
			// MOD-BY-LEETEN 2013/07/07-FROM:			pvpmapSparse = NULL;
			// ADD-BY-LEETEN 2013/07/23-BEGIN
			#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
			pvpmapSparse = NULL;
			#else	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		
			// ADD-BY-LEETEN 2013/07/23-END

			this->pcFullArrays = NULL;
			this->pcDecodingSparseArrays = NULL;
			this->pcEncodingSparseArrays = NULL;
			#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// ADD-BY-LEETEN 2013/07/23
			// MOD-BY-LEETEN 2013/07/07-END
		}

		// ADD-BY-LEETEN 11/12/2012-BEGIN
		virtual	// ADD-BY-LEETEN 01/02/2013
		~CSepDWTPool()
		{
			#if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:			#if	0	// MOD-BY-LEETEN 2013/07/07-FROM:
			if(this->pvpmapSparse)
			{	// ADD-BY-LEETEN 11/19/2012
				for(size_t e = 0; e < this->pvpmapSparse->size(); e++)
					if( (*this->pvpmapSparse)[e] )
					{
						(*this->pvpmapSparse)[e]->clear();
						delete (*this->pvpmapSparse)[e];
						(*this->pvpmapSparse)[e] = NULL;
					}
			// ADD-BY-LEETEN 11/19/2012-BEGIN
				delete [] this->pvpmapSparse;
				this->pvpmapSparse = NULL;
			}
			// ADD-BY-LEETEN 11/19/2012-END
			#else // #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#else	// MOD-BY-LEETEN 2013/07/07-TO:
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
			#endif	// #if	!WITH_DYNAMIC_ARRAY_ALLOCATION		// MOD-BY-LEETEN 2013/07/23-FROM:				#endif	// MOD-BY-LEETEN 2013/07/07-END
		}
		// ADD-BY-LEETEN 11/12/2012-END
	};
}
		
