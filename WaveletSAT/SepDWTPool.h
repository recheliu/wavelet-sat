#pragma once

// ADD-BY-LEETEN 11/12/2012-BEGIN
#define WITH_POINTER_TO_MAP	1
// ADD-BY-LEETEN 11/12/2012-END

#define WITH_SPARSE_AS_VECTOR	1	// ADD-BY-LEETEN 11/11/2012

#define	WITH_UNORDERED_MAP		1	// ADD-By-LEETEN 02/19/2013

#include	<algorithm>	// ADD-BY-LEETEN 03/17/2013

#if	!WITH_UNORDERED_MAP	// ADD-By-LEETEN 02/19/2013
#include <map>	
// ADD-By-LEETEN 02/19/2013-BEGIN
#else	// #if	!WITH_UNORDERED_MAP
#include <unordered_map>	
#endif	// #if	!WITH_UNORDERED_MAP
// ADD-By-LEETEN 02/19/2013-END
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

		vector< vector<WT> > vvFull;

		#if	!WITH_POINTER_TO_MAP	// ADD-BY-LEETE 11/12/2012
		vector< map<IT, ST> > vmapSparse;

		// ADD-BY-LEETEN 11/12/2012-BEGIN
		#else	// #if	!WITH_POINTER_TO_MAP
		#if	!WITH_UNORDERED_MAP	// ADD-By-LEETEN 02/19/2013
		vector< map<BT, WT>* >* pvpmapSparse;
		// ADD-By-LEETEN 02/19/2013-BEGIN
		#else	// #if	!WITH_UNORDERED_MAP	
		vector< unordered_map<BT, WT>* >* pvpmapSparse;
		#endif	// #if	!WITH_UNORDERED_MAP	
		// ADD-By-LEETEN 02/19/2013-END
		#endif	// #if	!WITH_POINTER_TO_MAP
		// ADD-BY-LEETEN 11/12/2012-END
		
		// ADD-BY-LEETEN 11/11/2012-BEGIN
		//! Max # that each coefficient is updated
		size_t uMaxCount;

		#if	WITH_SPARSE_AS_VECTOR
		//! Record how many time each coefficient has been updated
		vector<size_t> vuCounts;

		//! 
		vector< vector< pair<BT, WT> > > vvpairSparse;
		#endif	// #if	WITH_SPARSE_AS_VECTOR
		// ADD-BY-LEETEN 11/11/2012-END
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
		}
		// ADD-BY-LEETEN 01/27/2013-END

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
			if( !bIsSparse )
			{
				vvFull.resize(uNrOfBins);
				for(size_t b = 0; b < (size_t)uNrOfBins; b++)
					vvFull[b].resize(this->uSize);

			}
			else
			{	// ADD-BY-LEETEN 11/11/2012
				#if	!WITH_POINTER_TO_MAP	// ADD-BY-LEETEN 11/12/2012
				this->vmapSparse.resize(this->uSize);
				// ADD-BY-LEETEN 11/12/2012-BEGIN
				#else	// #if	!WITH_POINTER_TO_MAP	
				#if	!WITH_UNORDERED_MAP	// ADD-By-LEETEN 02/19/2013
				this->pvpmapSparse = new vector< map <BT, WT>* >;
				// ADD-By-LEETEN 02/19/2013-BEGIN
				#else	// #if	!WITH_UNORDERED_MAP	
				this->pvpmapSparse = new vector< unordered_map <BT, WT>* >;
				#endif	// #if	!WITH_UNORDERED_MAP	
				// ADD-By-LEETEN 02/19/2013-END
				this->pvpmapSparse->resize(this->uSize);
				#endif	// #if	!WITH_POINTER_TO_MAP	
				// ADD-BY-LEETEN 11/12/2012-END

			// ADD-BY-LEETEN 11/11/2012-BEGIN
				#if	WITH_SPARSE_AS_VECTOR
				this->vvpairSparse.resize(this->uSize);
				// DEL-BY-LEETEN 03/16/2013:	this->vuCounts.resize(this->uSize);
				#endif	// #if	WITH_SPARSE_AS_VECTOR
			}
			// ADD-BY-LEETEN 11/11/2012-END
			// MOD-BY-LEETEN 03/28/2013-FROM:			this->vuCounts.resize(this->uSize);	// ADD-BY-LEETEN 03/16/2013
			this->vuCounts.assign(this->uSize, 0);
			// MOD-BY-LEETEN 03/28/2013-END

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
				#if	!WITH_SPARSE_AS_VECTOR	// ADD-BY-LEETEN 11/11/2012
				#if	!WITH_POINTER_TO_MAP	// ADD-BY-LEETEN 11/12/2012
				for(typename vector< map<IT, ST> >::iterator 
					ivmapSparse = this->vmapSparse.begin();
					ivmapSparse != this->vmapSparse.end();
					ivmapSparse++)
				{
					map<IT, ST>& vmapBinSparse = *ivmapSparse;
					for(typename map<IT, ST>::iterator 
						pair = vmapBinSparse.begin();
						pair != vmapBinSparse.end();
						pair++)
						if( (ST)fabs(pair->second) >= Threshold )
							uCountInSparseArray++;
				}
				// ADD-BY-LEETEN 11/12/2012-BEGIN
				#else	// #if	!WITH_POINTER_TO_MAP	
				for(typename vector< map<IT, ST>* >::iterator 
					ivpmapSparse = this->pvpmapSparse->begin();
					ivpmapSparse != this->pvpmapSparse->end();
					ivpmapSparse++)
				{
					if(NULL == *ivpmapSparse)
						continue;

					map<IT, ST>& vmapBinSparse = *(*ivpmapSparse);
					for(typename map<IT, ST>::iterator 
						pair = vmapBinSparse.begin();
						pair != vmapBinSparse.end();
						pair++)
						if( (ST)fabs(pair->second) >= Threshold )
							uCountInSparseArray++;
				}
				#endif	// #if	!WITH_POINTER_TO_MAP	
				// ADD-BY-LEETEN 11/12/2012-END

				// ADD-BY-LEETEN 11/11/2012-BEGIN
				#else	// #if	!WITH_SPARSE_AS_VECTOR	
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
				#endif	// #if	!WITH_SPARSE_AS_VECTOR	
				// ADD-BY-LEETEN 11/11/2012-END
			}
		}



		void
		_AddEntryToSparseArray
		(
			#if	!WITH_UNORDERED_MAP	// ADD-By-LEETEN 02/19/2013
			map<BT, WT>& mapSparse,
			// ADD-By-LEETEN 02/19/2013-BEGIN
			#else	//	#if	!WITH_UNORDERED_MAP
			unordered_map<BT, WT>& mapSparse,
			#endif	//	#if	!WITH_UNORDERED_MAP	
			// ADD-By-LEETEN 02/19/2013-END
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
				Value = vvFull[Bin][uIndex];
			}
			else
			{
				// sparse
				vector< pair<BT, WT> >& vpair = this->vvpairSparse[uIndex];
				if( (size_t)usOffset < vpair.size() )
				{
					Bin = vpair[usOffset].first;
					Value = vpair[usOffset].second;
				}
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
				Value = vvFull[Bin][uIndex];
			}
			else
			{
				// sparse
				#if	!WITH_POINTER_TO_MAP	// ADD-BY-LEETEN 11/12/2012
				map<IT, ST>& mapCoef = vmapSparse[uIndex];

				// ADD-BY-LEETEN 11/12/2012-BEGIN
				#else	// #if	!WITH_POINTER_TO_MAP	
				if( NULL == (*pvpmapSparse)[uIndex] )
				{
					Value = WT(0);
					return;
				}

				#if	!WITH_UNORDERED_MAP	// ADD-By-LEETEN 02/19/2013
				map<BT, WT>& mapCoef = *(*pvpmapSparse)[uIndex];
				// ADD-By-LEETEN 02/19/2013-BEGIN
				#else	// #if	!WITH_UNORDERED_MAP	
				unordered_map<BT, WT>& mapCoef = *(*pvpmapSparse)[uIndex];
				#endif	// #if	!WITH_UNORDERED_MAP	
				// ADD-By-LEETEN 02/19/2013-END
				#endif	// #if	!WITH_POINTER_TO_MAP	
				// ADD-BY-LEETEN 11/12/2012-END
				#if	!WITH_UNORDERED_MAP	// ADD-By-LEETEN 02/19/2013
				typename map<BT, WT>::iterator ipair = mapCoef.find(Bin);
				// ADD-By-LEETEN 02/19/2013-BEGIN
				#else	// #if	!WITH_UNORDERED_MAP	
				typename unordered_map<BT, WT>::iterator ipair = mapCoef.find(Bin);
				#endif	// #if	!WITH_UNORDERED_MAP	
				// ADD-By-LEETEN 02/19/2013-END
				if( mapCoef.end() == ipair )
					Value = WT(0);
				else
					Value = ipair->second;
			}
		}

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
			#if	WITH_SPARSE_AS_VECTOR		
			// DEL-BY-LEETEN 03/16/2013:	if( bIsSparse )
				vuCounts[uIndex] += uCount;
			#endif	// #if	WITH_SPARSE_AS_VECTOR	

			if( Value )
			{
			// ADD-By-LEETEN 11/11/2012-END

			if( !bIsSparse )
			{
				// full
				vvFull[Bin][uIndex] += Value;
			}
			else
			{
				// sparse
				#if	!WITH_POINTER_TO_MAP	// ADD-BY-LEETEN 11/12/2012
				map<IT, ST>& mapCoef = vmapSparse[uIndex];

				// ADD-BY-LEETEN 11/12/2012-BEGIN
				#else	// #if	!WITH_POINTER_TO_MAP
				#if	!WITH_UNORDERED_MAP	// ADD-By-LEETEN 02/19/2013
				if(NULL == (*pvpmapSparse)[uIndex])
					(*pvpmapSparse)[uIndex] = new map<BT, WT>();
				map<BT, WT>& mapCoef = *(*pvpmapSparse)[uIndex];
				// ADD-By-LEETEN 02/19/2013-BEGIN
				#else	// #if	!WITH_UNORDERED_MAP	
				if(NULL == (*pvpmapSparse)[uIndex])
					(*pvpmapSparse)[uIndex] = new unordered_map<BT, WT>();
				unordered_map<BT, WT>& mapCoef = *(*pvpmapSparse)[uIndex];
				#endif	// #if	!WITH_UNORDERED_MAP	
				// ADD-By-LEETEN 02/19/2013-END
				#endif	// #if	!WITH_POINTER_TO_MAP	
				// ADD-BY-LEETEN 11/12/2012-END
				#if	!WITH_UNORDERED_MAP	// ADD-By-LEETEN 02/19/2013
				typename map<BT, WT>::iterator ipair = mapCoef.find(Bin);
				// ADD-By-LEETEN 02/19/2013-BEGIN
				#else	// #if	!WITH_UNORDERED_MAP	
				typename unordered_map<BT, WT>::iterator ipair = mapCoef.find(Bin);
				#endif	// #if	!WITH_UNORDERED_MAP	
				// ADD-By-LEETEN 02/19/2013-END
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
			}	// ADD-By-LEETEN 11/11/2012

			// ADD-BY-LEETEN 11/11/2012-BEGIN
			#if	WITH_SPARSE_AS_VECTOR	
			if( bIsSparse && vuCounts[uIndex] == uMaxCount )
			{
				#if	!WITH_POINTER_TO_MAP	// ADD-BY-LEETEN 11/12/2012
				const map<IT, ST>& vmapBinSparse = this->vmapSparse[uIndex];
				// ADD-BY-LEETEN 11/12/2012-BEGIN
				#else	// #if	!WITH_POINTER_TO_MAP	
				#if	!WITH_UNORDERED_MAP	// ADD-By-LEETEN 02/19/2013
				const map<BT, WT>& vmapBinSparse = *(*this->pvpmapSparse)[uIndex];
				// ADD-By-LEETEN 02/19/2013-BEGIN
				#else	// #if	!WITH_UNORDERED_MAP
				const unordered_map<BT, WT>& vmapBinSparse = *(*this->pvpmapSparse)[uIndex];
				#endif	// #if	!WITH_UNORDERED_MAP	
				// ADD-By-LEETEN 02/19/2013-END
				#endif	// #if	!WITH_POINTER_TO_MAP	
				// ADD-BY-LEETEN 11/12/2012-END

				vvpairSparse[uIndex].resize(vmapBinSparse.size());
				copy(vmapBinSparse.begin(), vmapBinSparse.end(), vvpairSparse[uIndex].begin());
				// ADD-By-LEETEN 02/19/2013-BEGIN
				#if		WITH_UNORDERED_MAP	
				sort(vvpairSparse[uIndex].begin(), vvpairSparse[uIndex].end());
				#endif	// #if	!WITH_UNORDERED_MAP	
				// ADD-By-LEETEN 02/19/2013-END
				/// now clear this map
				#if	!WITH_POINTER_TO_MAP	// ADD-BY-LEETEN 11/12/2012
				this->vmapSparse[uIndex].clear();

				// ADD-BY-LEETEN 11/12/2012-BEGIN
				#else	// #if	!WITH_POINTER_TO_MAP
				// _ShowMemoryUsage();
				(*this->pvpmapSparse)[uIndex]->clear();
				delete (*this->pvpmapSparse)[uIndex];
				(*this->pvpmapSparse)[uIndex] = NULL;
				// _ShowMemoryUsage();
				#endif	// #if	!WITH_POINTER_TO_MAP
				// ADD-BY-LEETEN 11/12/2012-END
			}	
			#endif	// #if	WITH_SPARSE_AS_VECTOR	
			// ADD-BY-LEETEN 11/11/2012-END
		}

		void
		_Finalize
		(
			WT WaveletWeight,
			void* _Reserved = NULL
		)
		{
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
				#if	!WITH_SPARSE_AS_VECTOR	// ADD-BY-LEETEN 11/11/2012
				#if	!WITH_POINTER_TO_MAP	// ADD-BY-LEETEN 11/12/2012
				for(typename vector< map<IT, ST> >::iterator 
					ivmapSparse = this->vmapSparse.begin();
					ivmapSparse != this->vmapSparse.end();
					ivmapSparse++)
				{
					map<IT, ST>& vmapBinSparse = *ivmapSparse;

				// ADD-BY-LEETEN 11/11/2012-BEGIN
				#else	// #if	!WITH_POINTER_TO_MAP
				for(typename vector< map<BT, WT>* >::iterator 
					ivpmapSparse = this->pvpmapSparse->begin();
					ivpmapSparse != this->pvpmapSparse->end();
					ivpmapSparse++)
				{
					map<BT, WT>& vmapBinSparse = *(*ivpmapSparse);
				#endif	// #if	!WITH_POINTER_TO_MAP	
				// ADD-BY-LEETEN 11/11/2012-END
					for(typename map<BT, WT>::iterator 
									pair = vmapBinSparse.begin();
						pair != vmapBinSparse.end();
						pair++)
						pair->second *= WaveletWeight;
				}
				// ADD-BY-LEETEN 11/11/2012-BEGIN
				#else	// #if	!WITH_SPARSE_AS_VECTOR	
				// check each element to make sure that no element is left in the map
				#if	!WITH_POINTER_TO_MAP	// ADD-BY-LEETEN 11/12/2012
				for(size_t e = 0; e < this->vmapSparse.size(); e++)
				{
					// ADD-By-LEETEN 11/11/2012-BEGIN
					// if the vector of pair is not empty, it means that the coefficients has been moved to here
					if(!vvpairSparse[e].empty())
						continue;
					// ADD-By-LEETEN 11/11/2012-END

					const map<IT, ST>& vmapBinSparse = this->vmapSparse[e];
					vvpairSparse[e].resize(vmapBinSparse.size());
					copy(vmapBinSparse.begin(), vmapBinSparse.end(), vvpairSparse[e].begin());

					/// now clear this map
					this->vmapSparse[e].clear();
				}
				// ADD-BY-LEETEN 11/12/2012-BEGIN
				#else	// #if	!WITH_POINTER_TO_MAP	
				for(size_t e = 0; e < this->pvpmapSparse->size(); e++)
				{
					// ADD-By-LEETEN 11/11/2012-BEGIN
					// if the vector of pair is not empty, it means that the coefficients has been moved to here
					if(NULL == (*this->pvpmapSparse)[e])
						continue;
					// ADD-By-LEETEN 11/11/2012-END

					#if		!WITH_UNORDERED_MAP	// ADD-By-LEETEN 02/19/2013
					const map<BT, WT>& vmapBinSparse = *(*this->pvpmapSparse)[e];
					// ADD-By-LEETEN 02/19/2013-BEGIN
					#else	// #if		!WITH_UNORDERED_MAP	
					const unordered_map<BT, WT>& vmapBinSparse = *(*this->pvpmapSparse)[e];
					#endif	// #if		!WITH_UNORDERED_MAP	
					// ADD-By-LEETEN 02/19/2013-END
					vvpairSparse[e].resize(vmapBinSparse.size());
					copy(vmapBinSparse.begin(), vmapBinSparse.end(), vvpairSparse[e].begin());
					// ADD-By-LEETEN 02/19/2013-BEGIN
					#if		WITH_UNORDERED_MAP	
					sort(vvpairSparse[e].begin(), vvpairSparse[e].end());
					#endif	//#if		WITH_UNORDERED_MAP	
					// ADD-By-LEETEN 02/19/2013-END
					/// now clear this map
					(*this->pvpmapSparse)[e]->clear();
					delete (*this->pvpmapSparse)[e];
					(*this->pvpmapSparse)[e] = NULL;
				}
				#endif	// #if	!WITH_POINTER_TO_MAP	
				// ADD-BY-LEETEN 11/12/2012-END

				// clear the vector to hold the sparse array and the counts
				#if	!WITH_POINTER_TO_MAP	// ADD-BY-LEETEN 11/12/2012
				this->vmapSparse.clear();

				// ADD-BY-LEETEN 11/12/2012-BEGIN
				#else	// #if	!WITH_POINTER_TO_MAP	
				this->pvpmapSparse->clear();
				delete this->pvpmapSparse;
				this->pvpmapSparse = NULL;
				#endif	// #if	!WITH_POINTER_TO_MAP	
				// ADD-BY-LEETEN 11/12/2012-END
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
				#endif	// #if	!WITH_SPARSE_AS_VECTOR	
				// ADD-BY-LEETEN 11/11/2012-END
			}
		}

		// ADD-BY-LEETEN 03/16/2013-BEGIN
		const 
		bool
		BIsEmpty(
			size_t uIndex,
			void* _Reserved = NULL
		) const
		{
			return (vuCounts[uIndex])?false:true;
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

		virtual
		const 
		vector< pair<BT, WT> >&
		VGetCoefSparse
		(
			size_t uIndex,
			void* _Reserved = NULL
		) const
		{
			if( !bIsSparse )
			{
				static vector< pair<BT, WT> > vpairCoefs;
				vpairCoefs.clear();

				for(size_t b = 0; b < this->vvFull.size(); b++)
				{
					WT Coef = vvFull[b][uIndex];
					if( Coef )
						vpairCoefs.push_back(pair<BT, WT>((BT)b, Coef));
				}
				return vpairCoefs;
			}
			else
				return this->vvpairSparse[uIndex];
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
				for(size_t b = 0; b < this->vvFull.size(); b++)
				{
					WT Coef = vvFull[b][uIndex];
					if( Coef )
						vpairCoefs.push_back(pair<BT, WT>((BT)b, Coef));
				}
			}
			else
			{
			// ADD-BY-LEETEN 11/12/2012-END

			// sparse
			#if	!WITH_SPARSE_AS_VECTOR	// ADD-BY-LEETEN 11/11/2012
			#if	!WITH_POINTER_TO_MAP	// ADD-BY-LEETEN 11/12/2012
			const map<IT, ST>& vmapBinSparse = this->vmapSparse[uIndex];

			// ADD-BY-LEETEN 11/12/2012-BEGIN
			#else	// #if	!WITH_POINTER_TO_MAP	
			const map<IT, ST>& vmapBinSparse = *(*this->pvpmapSparse)[uIndex];
			#endif	// #if	!WITH_POINTER_TO_MAP	
			// ADD-BY-LEETEN 11/12/2012-END

			vpairCoefs.clear();
			for(typename map<IT, ST>::const_iterator 
					ipair = vmapBinSparse.begin();
				ipair != vmapBinSparse.end();
				ipair++)
				vpairCoefs.push_back(pair<size_t, ST>((size_t)ipair->first, ipair->second));

			// ADD-BY-LEETEN 11/11/2012-BEGIN
			#else	// #if	!WITH_SPARSE_AS_VECTOR	
			const vector< pair<BT, WT> >& vpairSparse = this->vvpairSparse[uIndex];
			vpairCoefs.resize(vpairSparse.size());
			copy(vpairSparse.begin(), vpairSparse.end(), vpairCoefs.begin());
			#endif	// #if	!WITH_SPARSE_AS_VECTOR	
			// ADD-BY-LEETEN 11/11/2012-END
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
			// ADD-BY-LEETEN 11/12/2012-BEGIN
			#if	WITH_POINTER_TO_MAP	
			pvpmapSparse = NULL;
			#endif	// #if	WITH_POINTER_TO_MAP	
			// ADD-BY-LEETEN 11/12/2012-END
		}

		// ADD-BY-LEETEN 11/12/2012-BEGIN
		virtual	// ADD-BY-LEETEN 01/02/2013
		~CSepDWTPool()
		{
			#if	WITH_POINTER_TO_MAP	
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
			#endif	// #if	WITH_POINTER_TO_MAP
		}
		// ADD-BY-LEETEN 11/12/2012-END
	};
}
		
