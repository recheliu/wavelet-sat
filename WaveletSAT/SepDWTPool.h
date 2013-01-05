#pragma once

// ADD-BY-LEETEN 11/12/2012-BEGIN
#define WITH_POINTER_TO_MAP	1
// ADD-BY-LEETEN 11/12/2012-END

#define WITH_SPARSE_AS_VECTOR	1	// ADD-BY-LEETEN 11/11/2012

#include <map>	
#include <vector>
using namespace std;
#include <math.h>

#include "HeaderBase.h"	
#include "DWT.h"
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
	#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
	/*!
	ST: Type of the wavelet coefficients
	IT: Type of the bin index
	*/
	template<typename ST, typename IT = size_t>
	class CSepDWTPool
	#else	// MOD-BY-LEETEN 01/03/2013-TO:
	template<
		typename WT = typeWavelet,	//!< Type of the wavelet coefficientsd
		typename BT = typeBin		//!< Type of the bin
	>
	class CSepDWTPool
	#endif	// MOD-BY-LEETEN 01/03/2013-END
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
		vector< map<BT, WT>* >* pvpmapSparse;
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
		
		void
		_Set
		(
			// MOD-BY-LEETEN 01/03/2013-FROM:			const size_t uNrOfBins,
			const BT& uNrOfBins,
			// MOD-BY-LEETEN 01/03/2013-END
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
				// MOD-BY-LEETEN 01/03/2013-FROM:				for(size_t b = 0; b < uNrOfBins; b++)
				for(size_t b = 0; b < (size_t)uNrOfBins; b++)
				// MOD-BY-LEETEN 01/03/2013-END
					vvFull[b].resize(this->uSize);
			}
			else
			{	// ADD-BY-LEETEN 11/11/2012
				#if	!WITH_POINTER_TO_MAP	// ADD-BY-LEETEN 11/12/2012
				this->vmapSparse.resize(this->uSize);
				// ADD-BY-LEETEN 11/12/2012-BEGIN
				#else	// #if	!WITH_POINTER_TO_MAP	
				this->pvpmapSparse = new vector< map <BT, WT>* >;
				this->pvpmapSparse->resize(this->uSize);
				#endif	// #if	!WITH_POINTER_TO_MAP	
				// ADD-BY-LEETEN 11/12/2012-END

			// ADD-BY-LEETEN 11/11/2012-BEGIN
				#if	WITH_SPARSE_AS_VECTOR
				this->vvpairSparse.resize(this->uSize);
				this->vuCounts.resize(this->uSize);
				#endif	// #if	WITH_SPARSE_AS_VECTOR
			}
			// ADD-BY-LEETEN 11/11/2012-END

			this->uMaxCount = uMaxCount;	// ADD-BY-LEETEN 11/11/2012
		}

		void
		_GetArraySize
		(
			size_t& uCountInFullArray,
			size_t& uCountInSparseArray,
			// MOD-BY-LEETEN 01/03/2013-FROM:			ST Threshold,
			const WT& Threshold,
			// MOD-BY-LEETEN 01/03/2013-END
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
			#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
			map<IT, ST>& mapSparse,
			IT Bin,
			const ST& Value,
			#else	// MOD-BY-LEETEN 01/03/2013-TO:
			map<BT, WT>& mapSparse,
			const BT& Bin,
			const WT& Value,
			#endif	// MOD-BY-LEETEN 01/03/2013-END
			void* _Reserved = NULL
		)
		{
			if( Value )
			{
				mapSparse.insert(pair<BT, WT>(Bin, Value));

				#if	0	// DEL-BY-LEETEN 12/25/2012-BEGIN
				// estimate the current memory usage
				static size_t uCount;
				const size_t uMaxCount = 100000;
				if( 0 == uCount % uMaxCount )
				{
					LOG_VAR_TO_ERROR(uCount);
					_ShowMemoryUsage(true);
				}
				uCount++;
				#endif	// DEL-BY-LEETEN 12/25/2012-END
			}
		}

		// ADD-BY-LEETEN 12/29/2012-BEGIN
		void
		_GetAtOffset
		(
			// MOD-BY-LEETEN 01/03/2013-FROM:			const unsigned short usOffset,
			const BT& usOffset,
			// MOD-BY-LEETEN 01/03/2013-END
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

				map<BT, WT>& mapCoef = *(*pvpmapSparse)[uIndex];
				#endif	// #if	!WITH_POINTER_TO_MAP	
				// ADD-BY-LEETEN 11/12/2012-END
				typename map<BT, WT>::iterator ipair = mapCoef.find(Bin);
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
			void* _Reserved = NULL
		)
		{
			size_t uIndex = UConvertSubToIndex(vuSubs, vuLengths);
			// ADD-By-LEETEN 11/11/2012-BEGIN
			#if	WITH_SPARSE_AS_VECTOR		
			if( bIsSparse )
				vuCounts[uIndex]++;
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
				if(NULL == (*pvpmapSparse)[uIndex])
					(*pvpmapSparse)[uIndex] = new map<BT, WT>();
				map<BT, WT>& mapCoef = *(*pvpmapSparse)[uIndex];
				#endif	// #if	!WITH_POINTER_TO_MAP	
				// ADD-BY-LEETEN 11/12/2012-END
				typename map<BT, WT>::iterator ipair = mapCoef.find(Bin);
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
				const map<BT, WT>& vmapBinSparse = *(*this->pvpmapSparse)[uIndex];
				#endif	// #if	!WITH_POINTER_TO_MAP	
				// ADD-BY-LEETEN 11/12/2012-END

				vvpairSparse[uIndex].resize(vmapBinSparse.size());
				copy(vmapBinSparse.begin(), vmapBinSparse.end(), vvpairSparse[uIndex].begin());

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

					const map<BT, WT>& vmapBinSparse = *(*this->pvpmapSparse)[e];
					vvpairSparse[e].resize(vmapBinSparse.size());
					copy(vmapBinSparse.begin(), vmapBinSparse.end(), vvpairSparse[e].begin());

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

		void
		_GetCoefSparse
		(
			// MOD-BY-LEETEN 01/05/2013-FROM:			const vector<size_t> vuSubs,
			size_t uIndex,
			// MOD-BY-LEETEN 01/05/2013-END
			// MOD-BY-LEETEN 01/03/2013-FROM:			vector< pair<size_t, ST> >& vpairCoefs,
			vector< pair<BT, WT> >& vpairCoefs,
			// MOD-BY-LEETEN 01/03/2013-END
			void* _Reserved = NULL
		) const
		{
			// DEL-BY-LEETEN 01/05/2013:	size_t uIndex = UConvertSubToIndex(vuSubs, vuLengths);

			// ADD-BY-LEETEN 11/12/2012-BEGIN
			vpairCoefs.clear();
			if( !bIsSparse )
			{
				for(size_t b = 0; b < this->vvFull.size(); b++)
				{
					WT Coef = vvFull[b][uIndex];
					if( Coef )
						// MOD-BY-LEETEN 01/03/2013-FROM:						vpairCoefs.push_back(pair<size_t, ST>(b, Coef));
						vpairCoefs.push_back(pair<BT, WT>((BT)b, Coef));
						// MOD-BY-LEETEN 01/03/2013-END
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
			#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
			for(typename vector< pair<BT, WT> >::const_iterator 
							ipair = vpairSparse.begin();
				ipair != vpairSparse.end();
				ipair++)
				vpairCoefs.push_back(pair<size_t, ST>((size_t)ipair->first, ipair->second));
			#else	// MOD-BY-LEETEN 01/03/2013-TO:
			vpairCoefs.resize(vpairSparse.size());
			copy(vpairSparse.begin(), vpairSparse.end(), vpairCoefs.begin());
			#endif	// MOD-BY-LEETEN 01/03/2013-END
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
		
