#pragma once

#define WITH_SPARSE_AS_VECTOR	0	// ADD-BY-LEETEN 11/11/2012

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

	//! The base class of data pool of SepDWT coefficients. 
	/*!
	ST: Type of the wavelet coefficients
	IT: Type of the bin index
	*/
	template<typename ST, typename IT = size_t>
	class CSepDWTPool
	{
	protected:	
		// the pool of coefficients 
		//! Size of this coef. block. 
		vector<size_t> vuLengths;

		//! Size of this coef. block, which is equal to the product of vuLengths
		size_t uSize;

		//! A flag to indicate whether this coef. block is sparse.
		bool bIsSparse;

		vector< vector<ST> > vvFull;

		vector< map<IT, ST> > vmapSparse;
		
		// ADD-BY-LEETEN 11/11/2012-BEGIN
		//! Max # that each coefficient is updated
		size_t uMaxCount;

		#if	WITH_SPARSE_AS_VECTOR
		//! Record how many time each coefficient has been updated
		vector<size_t> vuCounts;

		//! 
		vector< vector< pair<IT, ST> > > vvpairSparse;
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
			const size_t uNrOfBins,
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
				for(size_t b = 0; b < uNrOfBins; b++)
					vvFull[b].resize(this->uSize);
			}
			else
			{	// ADD-BY-LEETEN 11/11/2012
				this->vmapSparse.resize(this->uSize);

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
			ST Threshold,
			void* _Reserved = NULL
		)
		{
			uCountInFullArray = 0;
			uCountInSparseArray = 0;
			if(!bIsSparse)
			{
				for(typename vector< vector<ST> >::iterator 
					ivvFull = this->vvFull.begin();
					ivvFull != this->vvFull.end();
					ivvFull++)
				{
					vector<ST>& vBinFull = *ivvFull;
					for(typename vector<ST>::iterator 
						ivFull = vBinFull.begin();
						ivFull != vBinFull.end();
						ivFull++)
						if( (ST)fabs(*ivFull) >= Threshold )
							uCountInFullArray++;
				}
			}
			else
			{
				#if	!WITH_SPARSE_AS_VECTOR	// ADD-BY-LEETEN 11/11/2012
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
				// ADD-BY-LEETEN 11/11/2012-BEGIN
				#else	// #if	!WITH_SPARSE_AS_VECTOR	
				for(typename vector< vector< pair<IT, ST> > >::iterator 
					ivvpairSparse = this->vvpairSparse.begin();
					ivvpairSparse != this->vvpairSparse.end();
					ivvpairSparse++)
				{
					vector< pair<IT, ST> >& vpairSparse = *ivvpairSparse;
					for(typename vector< pair<IT, ST> >::iterator 
						pair = vpairSparse.begin();
						pair != vpairSparse.end();
						pair++)
						if( (ST)fabs(pair->second) >= Threshold )
							uCountInSparseArray++;
				}
				#endif	// #if	!WITH_SPARSE_AS_VECTOR	
				// ADD-BY-LEETEN 11/11/2012-END
			}
		}



		void
		_AddEntryToSparseArray
		(
			map<IT, ST>& mapSparse,
			IT Bin,
			const ST& Value,
			void* _Reserved = NULL
		)
		{
			if( Value )
			{
				mapSparse.insert(pair<IT, ST>(Bin, Value));

				// estimate the current memory usage
				static size_t uCount;
				const size_t uMaxCount = 100000;
				if( 0 == uCount % uMaxCount )
				{
					LOG_VAR_TO_ERROR(uCount);
					_ShowMemoryUsage(true);
				}
				uCount++;
			}
		}

		void
		_GetAt
		(
			const IT Bin,
			const vector<size_t>& vuSubs,
			size_t& uIndex,
			ST& Value,
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
				map<IT, ST>& mapCoef = vmapSparse[uIndex];
				typename map<IT, ST>::iterator ipair = mapCoef.find(Bin);
				if( mapCoef.end() == ipair )
					Value = ST(0);
				else
					Value = ipair->second;
			}
		}

		#if	0	// DEL-BY-LEETEN 11/11/2012-BEGIN
		//! Set value to the location specified by the 1D index
		void
		_SetAt
		(
			const IT Bin,
			const vector<size_t>& vuSubs,
			const ST Value,
			void* _Reserved = NULL
		)
		{
			size_t uIndex = UConvertSubToIndex(vuSubs, vuLengths);
			if( !bIsSparse )
			{
				// full
				vvFull[Bin][uIndex] = Value;
			}
			else
			{
				// sparse
				map<IT, ST>& mapCoef = vmapSparse[uIndex];
				typename map<IT, ST>::iterator ipair = mapCoef.find(Bin);
				if( mapCoef.end() == ipair )
					_AddEntryToSparseArray
					(
						mapCoef,
						Bin,
						Value
					);
				else
					ipair->second = Value;
			}
		}
		#endif	// DEL-BY-LEETEN 11/11/2012-END

		//! Add value to the location specified by the 1D index
		void
		_AddAt
		(
			const IT Bin,
			const vector<size_t>& vuSubs,
			const ST Value,
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
				map<IT, ST>& mapCoef = vmapSparse[uIndex];
				typename map<IT, ST>::iterator ipair = mapCoef.find(Bin);
				if( mapCoef.end() == ipair )
					_AddEntryToSparseArray
					(
						mapCoef,
						Bin,
						Value
					);
					else
						ipair->second += Value;

				#if	0	// DEL-BY-LEETEN 11/11/2012-BEGIN
				// ADD-BY-LEETEN 11/11/2012-BEGIN
				#if	WITH_SPARSE_AS_VECTOR		
				vuCounts[uIndex]++;
				#endif	// #if	WITH_SPARSE_AS_VECTOR	
				// ADD-BY-LEETEN 11/11/2012-END
				#endif		// DEL-BY-LEETEN 11/11/2012-END
			}
			}	// ADD-By-LEETEN 11/11/2012

			// ADD-BY-LEETEN 11/11/2012-BEGIN
			#if	WITH_SPARSE_AS_VECTOR	
			if( bIsSparse && vuCounts[uIndex] == uMaxCount )
			{
				const map<IT, ST>& vmapBinSparse = this->vmapSparse[uIndex];
				#if	0	// MOD-BY-LEETEN 11/11/2012-FROM:
				vvpairSparse[uIndex].clear();
				for(typename map<IT, ST>::const_iterator 
					ipair = vmapBinSparse.begin();
					ipair != vmapBinSparse.end();
					ipair++)
					vvpairSparse[uIndex].push_back(pair<IT, ST>(ipair->first, ipair->second));
				#else		// MOD-BY-LEETEN 11/11/2012-TO:
				vvpairSparse[uIndex].resize(vmapBinSparse.size());
				copy(vmapBinSparse.begin(), vmapBinSparse.end(), vvpairSparse[uIndex].begin());
				#endif		// MOD-BY-LEETEN 11/11/2012-END
				/// now clear this map
				this->vmapSparse[uIndex].clear();
			}	
			#endif	// #if	WITH_SPARSE_AS_VECTOR	
			// ADD-BY-LEETEN 11/11/2012-END
		}

		void
		// MOD-BY-LEETEN 11/11/2012-FROM:	_Weight
		_Finalize
		// MOD-BY-LEETEN 11/11/2012-END
		(
			ST WaveletWeight,
			void* _Reserved = NULL
		)
		{
			if(!bIsSparse)
			{
				for(typename vector< vector<ST> >::iterator 
					ivvFull = this->vvFull.begin();
					ivvFull != this->vvFull.end();
					ivvFull++)
				{
					vector<ST>& vBinFull = *ivvFull;
					for(typename vector<ST>::iterator 
						ivFull = vBinFull.begin();
						ivFull != vBinFull.end();
						ivFull++)
					{
						ST Coef = *ivFull;
						if(Coef)
							*ivFull *= WaveletWeight;
					}
				}
			}
			else
			{
				#if	!WITH_SPARSE_AS_VECTOR	// ADD-BY-LEETEN 11/11/2012
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
						pair->second *= WaveletWeight;
				}
				// ADD-BY-LEETEN 11/11/2012-BEGIN
				#else	// #if	!WITH_SPARSE_AS_VECTOR	
				// check each element to make sure that no element is left in the map
				for(size_t e = 0; e < this->vmapSparse.size(); e++)
				{
					// ADD-By-LEETEN 11/11/2012-BEGIN
					// if the vector of pair is not empty, it means that the coefficients has been moved to here
					if(!vvpairSparse[e].empty())
						continue;
					// ADD-By-LEETEN 11/11/2012-END

					const map<IT, ST>& vmapBinSparse = this->vmapSparse[e];
					#if	0	// MOD-BY-LEETEN 11/11/2012-FROM:
					vvpairSparse[e].clear();
					for(typename map<IT, ST>::const_iterator 
						ipair = vmapBinSparse.begin();
						ipair != vmapBinSparse.end();
						ipair++)
					{
						vvpairSparse[e].push_back(pair<IT, ST>(ipair->first, ipair->second));
					}
					#else		// MOD-BY-LEETEN 11/11/2012-TO:
					vvpairSparse[e].resize(vmapBinSparse.size());
					copy(vmapBinSparse.begin(), vmapBinSparse.end(), vvpairSparse[e].begin());
					#endif		// MOD-BY-LEETEN 11/11/2012-END

					/// now clear this map
					this->vmapSparse[e].clear();
				}

				// clear the vector to hold the sparse array and the counts
				this->vmapSparse.clear();
				this->vuCounts.clear();

				for(typename vector< vector< pair<IT, ST> > >::iterator 
					ivvpairSparse = this->vvpairSparse.begin();
					ivvpairSparse != this->vvpairSparse.end();
					ivvpairSparse++)
				{
					vector< pair<IT, ST> >& vpairSparse = *ivvpairSparse;
					for(typename vector< pair<IT, ST> >::iterator 
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
			const vector<size_t> vuSubs,
			vector< pair<size_t, ST> >& vpairCoefs,
			void* _Reserved = NULL
		) const
		{
			size_t uIndex = UConvertSubToIndex(vuSubs, vuLengths);

			// sparse
			#if	!WITH_SPARSE_AS_VECTOR	// ADD-BY-LEETEN 11/11/2012
			const map<IT, ST>& vmapBinSparse = this->vmapSparse[uIndex];
			vpairCoefs.clear();
			for(typename map<IT, ST>::const_iterator 
					ipair = vmapBinSparse.begin();
				ipair != vmapBinSparse.end();
				ipair++)
				vpairCoefs.push_back(pair<size_t, ST>((size_t)ipair->first, ipair->second));

			// ADD-BY-LEETEN 11/11/2012-BEGIN
			#else	// #if	!WITH_SPARSE_AS_VECTOR	
			const vector< pair<IT, ST> >& vpairSparse = this->vvpairSparse[uIndex];
			#if	0	// MOD-BY-LEETEN 11/11/2012-FROM:
			for(typename vector< pair<IT, ST>>::const_iterator 
			#else		// MOD-BY-LEETEN 11/11/2012-TO:
			for(typename vector< pair<IT, ST> >::const_iterator 
			#endif		// MOD-BY-LEETEN 11/11/2012-END
					ipair = vpairSparse.begin();
				ipair != vpairSparse.end();
				ipair++)
				vpairCoefs.push_back(pair<size_t, ST>((size_t)ipair->first, ipair->second));
			#endif	// #if	!WITH_SPARSE_AS_VECTOR	
			// ADD-BY-LEETEN 11/11/2012-END
		}

		CSepDWTPool()
		{
		}
	};
}
		