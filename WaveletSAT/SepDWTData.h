#pragma once

//! Decide whether a table of the product of wavelet bais is precomputed.
/*!
This table is used in the decompression stage.
*/
#define	WITH_PRECOMPUTED_WAVELET_BASIS	0	

//! Decide whether a table of the product of wavelet sum is precomputed.
/*!
This table is used in the compression stage.
*/
#define	WITH_PRECOMPUTED_WAVELET_SUMS	0	

//! Decide whether a table that immediately map each updating coefficients and its dimension to the corresponding 1D index in the wavelet table per dimension.
#define WITH_COEF_DIM_2_WAVELET		0
// ADD-BY-LEETEN 10/21/2012-END

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
	template<typename ST>
	class CSepDWTData
	{
protected:	
		//! The maximun number of wavelet levels per dim, which are specified by the user.
		/*!
		m[0], ... m[d], ... m[D - 1]
		*/
		vector<size_t> vuDimMaxLevels;
		
		//! #levels per dim
		/*!
		l[0] = log2(n[0]) + 1, ... l[d] = log2(n[d]) + 1, ... l[D - 1] = log2(n[D - 1]) + 1
		*/
		vector<size_t> vuDimLevels;

		vector<size_t> vuDimLengths;
		vector<size_t> vuFullArrayDimLengths;
		vector<ST> vFullArray;
		map<size_t, ST> mapSparseArray;
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0700,
			PARAMETER_END
		};

		void
		_Set
		(
			const vector<size_t>& vuDimLengths,
			const size_t uMaxFullArraySize,
			void* _Reserved = NULL
		)
		{
			size_t uNrOfDims = vuDimLengths.size();
			size_t uMaxLevel = (size_t)floor((log((double)uMaxFullArraySize) / M_LN2) / (double)uNrOfDims);
			size_t uFullArraySize = (size_t)1<<(uMaxLevel * uNrOfDims);
			size_t uFullArrayDimLength = (size_t)1<<uMaxLevel;
			vuFullArrayDimLengths.resize(uNrOfDims);
			for(size_t d = 0; d < uNrOfDims; d++)
				vuFullArrayDimLengths[d] = uFullArrayDimLength;

			this->vuDimLengths.resize(uNrOfDims);
			for(size_t d = 0; d < uNrOfDims; d++)
				this->vuDimLengths[d] = vuDimLengths[d];
	
			vFullArray.resize(uFullArraySize);
		}

		//! Given a position, decide whether it is within the full array
		bool BIsInFullArray
		(
			const vector<size_t>& vuPos,
			void* _Reserved = NULL
		)
		{
			bool bIsInFullArray = true;
			for(size_t d = 0; d < vuPos.size(); d++)
				if( vuPos[d] >= vuFullArrayDimLengths[d] )
				{
					bIsInFullArray = false;
					break;
				}

			return bIsInFullArray;

		}

		void
		_AddEntryToSparseArray
		(
			const size_t uIndex,
			const ST& Value,
			void* _Reserved = NULL
		)
		{
			if( Value )
			{
				mapSparseArray.insert(pair<size_t, ST>(uIndex, Value));

				// estimate the current memory usage
				static size_t uCount;
				const size_t uMaxCount = 100000;
				if( 0 == uCount % uMaxCount )
				{
					LOG_VAR(uCount);
					_ShowMemoryUsage();
				}
				uCount++;
			}
		}

		//! Get the value at a position
		void
		_GetAtPos
		(
			const vector<size_t>& vuPos,
			ST& Value,
			void* _Reserved = NULL
		)
		{
			bool bIsInFullArray = BIsInFullArray(vuPos);
			if(bIsInFullArray)
			{
				size_t uIndexInFullArray = UConvetSubToIndex(vuPos, vuFullArrayDimLengths);
				Value = vFullArray[uIndexInFullArray];
			}
			else
			{
				size_t uIndex = UConvetSubToIndex(vuPos, vuDimLengths);
				map<size_t, double>::iterator ipair = mapSparseArray.find(uIndex);
				if(mapSparseArray.end() == ipair )
					Value = 0;
				else
					Value = ipair->second;
			}
		}

		//! Set the value at a position
		void
		_SetAtPos
		(
			const vector<size_t>& vuPos,
			const ST Value,
			void* _Reserved = NULL
		)
		{
			bool bIsInFullArray = BIsInFullArray(vuPos);
			if(bIsInFullArray)
			{
				size_t uIndexInFullArray = UConvetSubToIndex(vuPos, vuFullArrayDimLengths);
				vFullArray[uIndexInFullArray] = Value;
			}
			else
			{
				size_t uIndex = UConvetSubToIndex(vuPos, vuDimLengths);
				map<size_t, double>::iterator ipair = mapSparseArray.find(uIndex);
				if( mapSparseArray.end() == ipair )
					_AddEntryToSparseArray(uIndex, Value);
				else
					ipair->second = Value;
			}
		}

		//! Set the value at a position
		void
		_AddAtPos
		(
			const vector<size_t>& vuPos,
			const ST Value,
			void* _Reserved = NULL
		)
		{
			bool bIsInFullArray = BIsInFullArray(vuPos);
			if(bIsInFullArray)
			{
				size_t uIndexInFullArray = UConvetSubToIndex(vuPos, vuFullArrayDimLengths);
				vFullArray[uIndexInFullArray] += Value;
			}
			else
			{
				size_t uIndex = UConvetSubToIndex(vuPos, vuDimLengths);
				map<size_t, double>::iterator ipair = mapSparseArray.find(uIndex);
				if( mapSparseArray.end() == ipair )
					_AddEntryToSparseArray(uIndex, Value);
				else
					ipair->second += Value;
			}
		}

		//! Get the value to the location specified by the 1D index
		void
		_GetAtIndex
		(
			const size_t uIndex,
			ST& Value,
			void* _Reserved = NULL
		)
		{
#if	0 // TMP-MOD
			vector<size_t> vuPos;
			_ConvertIndexToSub(uIndex, vuPos, vuDimLengths);
			_GetAtPos(vuPos, Value);
#else
			Value = this->vFullArray[uIndex];
#endif
		}

		//! Set value to the location specified by the 1D index
		void
		_SetAtIndex
		(
			const size_t uIndex,
			const ST Value,
			void* _Reserved = NULL
		)
		{
#if	0	// TMP-MOD
			vector<size_t> vuPos;
			_ConvertIndexToSub(uIndex, vuPos, vuDimLengths);
			_SetAtPos(vuPos, Value);
#else
			vFullArray[uIndex] = Value;
#endif
		}

		//! Add value to the location specified by the 1D index
		void
		_AddAtIndex
		(
			const size_t uIndex,
			const ST Value,
			void* _Reserved = NULL
		)
		{
#if	0 // TMP-MOD
			vector<size_t> vuPos;
			_ConvertIndexToSub(uIndex, vuPos, vuDimLengths);
			_AddAtPos(vuPos, Value);
#else
			this->vFullArray[uIndex] += Value;
#endif
		}

		void
		_Finalize
		(
			ST WaveletDenomiator,
			void* _Reserved = NULL
		)
		{
			vector<size_t> vuSub;	// ADD-BY-LEETEN 10/06/2012
			for(size_t f = 0; f < this->vFullArray.size(); f++)
			{
				ST dCoef = this->vFullArray[f];
				if( dCoef )
				{
					ST Wavelet = (ST)+1.0;
					_ConvetIndexToSub(f, vuSub, vuDimLengths);
					for(size_t d = 0; d < vuSub.size(); d++)
					{
						size_t uSub = vuSub[d];
						if( uSub >= 1 )
						{
							size_t uLevel = (size_t)ceil(log( (double)(uSub + 1) ) / M_LN2 );
							Wavelet *= (ST)sqrt((double)(1 << (uLevel - 1) ));
						}
					}
					this->vFullArray[f] *= Wavelet / WaveletDenomiator;
				}
			}

			for(map<size_t, double>::iterator 
				ipairCoef = this->mapSparseArray.begin();
				ipairCoef != this->mapSparseArray.end();
				ipairCoef++)
			{
				double dCoef = ipairCoef->second;
				if( dCoef )
				{
					ST Wavelet = 1.0;

					_ConvetIndexToSub(ipairCoef->first, vuSub, vuDimLengths);
					for(size_t d = 0; d < vuSub.size(); d++)
					{
						size_t uSub = vuSub[d];
						if( uSub >= 1 )
						{
							size_t uLevel = (size_t)ceil(log( (double)(uSub + 1) ) / M_LN2 );
							Wavelet *= (ST)sqrt((double)(1 << (uLevel - 1) ));	
						}
					}
					ipairCoef->second *= Wavelet / WaveletDenomiator;
				}
			}
		}

		//! Add value to the location specified by the 1D index
		void
		_GetNrOfNonZeroCoefs
		(
			size_t& uCount,
			ST Threshold,
			void* _Reserved = NULL
		)
		{
			uCount = 0;
			for(size_t w = 0; w < this->vFullArray.size(); w++)
			{
				double dCoef = this->vFullArray[w];
				if( fabs(dCoef) > Threshold )
					uCount++;
			}

			for(map<size_t, ST>::iterator
				ipairCoef = this->mapSparseArray.begin();
				ipairCoef != this->mapSparseArray.end();
				ipairCoef++)
			{
				double dCoef = ipairCoef->second;
				if( fabs(dCoef) > Threshold )
					uCount++;
			}
		}

		CSepDWTData()
		{
		}
	};
}
