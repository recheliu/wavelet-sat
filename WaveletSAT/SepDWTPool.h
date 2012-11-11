#pragma once

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
				this->vmapSparse.resize(this->uSize);
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
			}
		}

		void
		_Weight
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
			const map<IT, ST>& vmapBinSparse = this->vmapSparse[uIndex];
			vpairCoefs.clear();
			for(typename map<IT, ST>::const_iterator 
					ipair = vmapBinSparse.begin();
				ipair != vmapBinSparse.end();
				ipair++)
				vpairCoefs.push_back(pair<size_t, ST>((size_t)ipair->first, ipair->second));
		}

		CSepDWTPool()
		{
		}
	};
}
		