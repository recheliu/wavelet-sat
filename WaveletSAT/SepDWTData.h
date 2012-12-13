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

// ADD-BY-LEETEN 10/21/2012-END

// ADD-BY-LEETEN 10/31/2012-BEGIN
//! Whether the full array and sparse array are divided based on 1D index to the multi-dimensional subscripts
#define WITH_1D_DIVISION		1
// ADD-BY-LEETEN 10/31/2012-END

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

	//! SepDWT coefficients per basis and per bin.
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
		#if	!WITH_1D_DIVISION		// ADD-BY-LEETEN 10/31/2012
		vector<size_t> vuFullArrayDimLengths;
		#endif	// #if	!WITH_1D_DIVISION	// ADD-BY-LEETEN 10/31/2012
		vector<ST> vFullArray;
		map<size_t, ST> mapSparseArray;

		bool bIsFullArrayOnly;	// ADD-BY-LEETEN 10/30/2012
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

		/*
		uMaxFullArraySize: in bytes
		*/
		void
		_Set
		(
			const vector<size_t>& vuDimLengths,
			const size_t uMaxFullArraySize,
			void* _Reserved = NULL
		)
		{
			size_t uNrOfDims = vuDimLengths.size();
			#if	!WITH_1D_DIVISION	// ADD-BY-LEETEN 10/31/2012
			vuDimLevels.resize(uNrOfDims);
			size_t uLevelsProduct = 1;	// product of levels from all dim
			for(size_t d = 0; d < uNrOfDims; d++)
			{
				size_t uNrOfDimLevels = (size_t)ceil((log((double)vuDimLengths[d]) / M_LN2)) + 1;
				vuDimLevels[d] = uNrOfDimLevels;
				uLevelsProduct *= uNrOfDimLevels;
			}

			vuFullArrayDimLengths.resize(uNrOfDims);
			size_t uFullArraySize = 0;
			if(uMaxFullArraySize)
			{
				vector<size_t> vuOptimalDimLevel;
				vuOptimalDimLevel.resize(uNrOfDims);

				size_t uDiff = uMaxFullArraySize;
				double dCurrentAspectRatio = -1.0; // ADD-BY-LEETEN 10/31/2012
				for(size_t l = 0; l < uLevelsProduct; l++)
				{
					vector<size_t> vuLevel;
					_ConvertIndexToSub(l, vuLevel, vuDimLevels);
					size_t uSize = 1;	// sizeof(ST);

					// ADD-BY-LEETEN 10/31/2012-BEGIN
					double dMaxLength = -HUGE_VAL;
					double dMinLength = +HUGE_VAL;
					// ADD-BY-LEETEN 10/31/2012-END

					for(size_t d = 0; d < uNrOfDims; d++)
					{
					  size_t uLength = (size_t)(1 << vuLevel[d]);
					  uSize *= uLength;
					  dMaxLength = max(dMaxLength, (double)uLength);
					  dMinLength = min(dMinLength, (double)uLength);
					}
					double dAspectRatio = dMaxLength / dMinLength;

					if( uSize <= uMaxFullArraySize )
					{
						size_t uNewDiff = uMaxFullArraySize - uSize;
						if( uNewDiff <= uDiff )
						  if( uNewDiff < uDiff ||
						      (dCurrentAspectRatio < 0.0 || dAspectRatio <= dCurrentAspectRatio ) )
						    {
							vuOptimalDimLevel = vuLevel;
							uDiff = uNewDiff;
							dCurrentAspectRatio = dAspectRatio;
						    }
					}
				}

				uFullArraySize = 1;
				for(size_t d = 0; d < uNrOfDims; d++)
				{
					vuFullArrayDimLengths[d] = (size_t)(1 << vuOptimalDimLevel[d]);
					uFullArraySize *= vuFullArrayDimLengths[d];
				}
			}

			this->vuDimLengths.resize(uNrOfDims);
			bIsFullArrayOnly = true;
			for(size_t d = 0; d < uNrOfDims; d++)
			{
				this->vuDimLengths[d] = vuDimLengths[d];
				if( this->vuFullArrayDimLengths[d] < this->vuDimLengths[d] )
					bIsFullArrayOnly = false;
			}
			// ADD-BY-LEETEN 10/31/2012-BEGIN
			#else	// #if	!WITH_1D_DIVISION	
			vuDimLevels.resize(uNrOfDims);
			for(size_t d = 0; d < uNrOfDims; d++)
			{
				size_t uNrOfDimLevels = (size_t)ceil((log((double)vuDimLengths[d]) / M_LN2)) + 1;
				vuDimLevels[d] = uNrOfDimLevels;
			}

			this->vuDimLengths.resize(uNrOfDims);
			size_t uDataSize = 1;
			for(size_t d = 0; d < uNrOfDims; d++)
			{
				this->vuDimLengths[d] = vuDimLengths[d];
				uDataSize *= vuDimLengths[d];
			}
			bIsFullArrayOnly = ( uDataSize < uMaxFullArraySize )?true:false;
			size_t uFullArraySize = min(uDataSize, uMaxFullArraySize);
			#endif	// #if	!WITH_1D_DIVISION	
			// ADD-BY-LEETEN 10/31/2012-END

			vFullArray.resize(uFullArraySize);
		}

#if !WITH_1D_DIVISION	       // ADD-BY-LEETEN 10/31/2012
		//! Given a position, decide whether it is within the full array
		bool BIsInFullArray
		(
			const vector<size_t>& vuPos,
			void* _Reserved = NULL
		)
		{
			bool bIsInFullArray = true;
			if( !bIsFullArrayOnly )	// ADD-BY-LEETEN 10/30/2012
			for(size_t d = 0; d < vuPos.size(); d++)
				if( vuPos[d] >= vuFullArrayDimLengths[d] )
				{
					bIsInFullArray = false;
					break;
				}

			return bIsInFullArray;

		}
#endif // #if !WITH_1D_DIVISION	       // ADD-BY-LEETEN 10/31/2012

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
					LOG_VAR_TO_ERROR(uCount);
					_ShowMemoryUsage(true);
				}
				uCount++;
			}
		}

		#if	!WITH_1D_DIVISION	// ADD-BY-LEETEN 10/31/2012
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
				size_t uIndexInFullArray = UConvertSubToIndex(vuPos, vuFullArrayDimLengths);
				Value = vFullArray[uIndexInFullArray];
			}
			else
			{
				size_t uIndex = UConvertSubToIndex(vuPos, vuDimLengths);
				typename map<size_t, ST>::iterator ipair = mapSparseArray.find(uIndex);
				if(mapSparseArray.end() == ipair )
					Value = ST(0);
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
				size_t uIndexInFullArray = UConvertSubToIndex(vuPos, vuFullArrayDimLengths);
				vFullArray[uIndexInFullArray] = Value;
			}
			else
			{
				size_t uIndex = UConvertSubToIndex(vuPos, vuDimLengths);
				typename map<size_t, ST>::iterator ipair = mapSparseArray.find(uIndex);
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
				size_t uIndexInFullArray = UConvertSubToIndex(vuPos, vuFullArrayDimLengths);
				vFullArray[uIndexInFullArray] += Value;
			}
			else
			{
				size_t uIndex = UConvertSubToIndex(vuPos, vuDimLengths);
				typename map<size_t, ST>::iterator ipair = mapSparseArray.find(uIndex);
				if( mapSparseArray.end() == ipair )
					_AddEntryToSparseArray(uIndex, Value);
				else
					ipair->second += Value;
			}
		}

		// ADD-BY-LEETEN 10/31/2012-BEGIN
		#else	// #if	!WITH_1D_DIVISION	
		//! Get the value to the location specified by the 1D index
		void
		_GetAtIndex
		(
			const size_t uIndex,
			ST& Value,
			void* _Reserved = NULL
		)
		{
			if( uIndex < vFullArray.size() )
				Value = vFullArray[uIndex];
			else
			{
				typename map<size_t, ST>::iterator ipair = mapSparseArray.find(uIndex);
				if(mapSparseArray.end() == ipair )
					Value = ST(0);
				else
					Value = ipair->second;
			}
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
			if( uIndex < vFullArray.size() )
				vFullArray[uIndex] = Value;
			else
			{
				typename map<size_t, ST>::iterator ipair = mapSparseArray.find(uIndex);
				if( mapSparseArray.end() == ipair )
					_AddEntryToSparseArray(uIndex, Value);
				else
					ipair->second = Value;
			}
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
			if( uIndex < vFullArray.size() )
				vFullArray[uIndex] += Value;
			else
			{
				typename map<size_t, ST>::iterator ipair = mapSparseArray.find(uIndex);
				if( mapSparseArray.end() == ipair )
					_AddEntryToSparseArray(uIndex, Value);
				else
					ipair->second += Value;
			}
		}
		#endif	// #if	!WITH_1D_DIVISION	
		// ADD-BY-LEETEN 10/31/2012-END

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
					#if !WITH_1D_DIVISION	// ADD-BY-LEETEN 10/31/2012
					_ConvertIndexToSub(f, vuSub, this->vuFullArrayDimLengths);
					// ADD-BY-LEETEN 10/31/2012-BEGIN
					#else	// #if !WITH_1D_DIVISION	
					_ConvertIndexToSub(f, vuSub, this->vuDimLengths);
					#endif	// #if !WITH_1D_DIVISION	
					// ADD-BY-LEETEN 10/31/2012-END
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

			for(typename map<size_t, ST>::iterator 
				ipairCoef = this->mapSparseArray.begin();
				ipairCoef != this->mapSparseArray.end();
				ipairCoef++)
			{
				ST dCoef = ipairCoef->second;
				if( dCoef )
				{
					ST Wavelet = 1.0;

					_ConvertIndexToSub(ipairCoef->first, vuSub, vuDimLengths);
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
		_GetArraySize
		(
			size_t& uCountInFullArray,
			size_t& uCountInSparseArray,
			ST Threshold,
			void* _Reserved = NULL
		)
		{
			uCountInFullArray = 0;
			for(size_t w = 0; w < this->vFullArray.size(); w++)
			{
				double dCoef = (double)this->vFullArray[w];
				if( fabs(dCoef) > Threshold )
					uCountInFullArray++;
			}

			uCountInSparseArray = 0;
			for(typename map<size_t, ST>::iterator
				ipairCoef = this->mapSparseArray.begin();
				ipairCoef != this->mapSparseArray.end();
				ipairCoef++)
			{
				double dCoef = (double)ipairCoef->second;
				if( fabs(dCoef) > Threshold )
					uCountInSparseArray++;
			}
		}

		CSepDWTData()
		{
		}
	};
}
