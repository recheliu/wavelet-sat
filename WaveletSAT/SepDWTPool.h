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

		vector< vector<WT> > vvFull;

		vector< unordered_map<BT, WT>* >* pvpmapSparse;
		
		// ADD-BY-LEETEN 11/11/2012-BEGIN
		//! Max # that each coefficient is updated
		size_t uMaxCount;

		//! Record how many time each coefficient has been updated
		vector<size_t> vuCounts;

		//! 
		vector< vector< pair<BT, WT> > > vvpairSparse;
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
				this->pvpmapSparse = new vector< unordered_map <BT, WT>* >;
				this->pvpmapSparse->resize(this->uSize);

			// ADD-BY-LEETEN 11/11/2012-BEGIN
				this->vvpairSparse.resize(this->uSize);
			}
			// ADD-BY-LEETEN 11/11/2012-END
			this->vuCounts.assign(this->uSize, 0);

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
		}



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
				vuCounts[uIndex] += uCount;

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
			}	// ADD-By-LEETEN 11/11/2012

			if( bIsSparse && vuCounts[uIndex] == uMaxCount )
			{
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
			const vector< pair<BT, WT> >& vpairSparse = this->vvpairSparse[uIndex];
			vpairCoefs.resize(vpairSparse.size());
			copy(vpairSparse.begin(), vpairSparse.end(), vpairCoefs.begin());
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
			pvpmapSparse = NULL;
		}

		// ADD-BY-LEETEN 11/12/2012-BEGIN
		virtual	// ADD-BY-LEETEN 01/02/2013
		~CSepDWTPool()
		{
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
		}
		// ADD-BY-LEETEN 11/12/2012-END
	};
}
		
