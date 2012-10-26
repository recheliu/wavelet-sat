#pragma once

#include <iostream>
#include <vector>
using namespace std;
#include <math.h>

#include "liblog.h"

#include "EncoderBase.h"

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

	//! The class that create SAT in core.
	/*!
	*/
	template<typename DT, typename ST>
	class CSATEncoder:
		public CHeaderBase,
		public CEncoderBase<DT, ST>
	{
protected:	
		//! The storage to store the original data.
		vector< vector<ST> > vvBinSATs;
		
		////////////////////////////////////////////////////////////////////
		/*
		The protected interface. 
		*/
		//! Update the specified bin.
		virtual
		void 
		_UpdateBin
		(
			const vector<size_t>& vuPos, 
			const DT& value,
			size_t uBin, 
			const ST& weight,
			void *_Reserved = NULL
		)
		{
			size_t uIndex = UConvetSubToIndex(vuPos, vuDimLengths);
			this->vvBinSATs[uBin][uIndex] += weight;
		}

public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/
		enum EParameter
		{
			PARAMETER_BEGIN = 0x0300,
			PARAMETER_END
		};

		//! Finalize the computation of SAT
		/*!
		*/
		virtual 
		void 
		_Finalize
		(
			void *_Reserved = NULL
		)
		{
			vector<size_t> vuScanLineBase;	// the base coordinate of the scane line
			vuScanLineBase.resize(this->UGetNrOfDims());

			vector<size_t> vuOtherDimLengths;	// the base coordinate of the scane line
			vuOtherDimLengths.resize(this->UGetNrOfDims());
			
			for(size_t uOffset = 1, d = 0; d < this->UGetNrOfDims(); uOffset *= vuDimLengths[d], d++)
			{
				size_t uNrOfScanLines = this->uDataSize / this->vuDimLengths[d];
				vuOtherDimLengths = vuDimLengths;
				vuOtherDimLengths[d] = 1;

				for(size_t i = 0; i < uNrOfScanLines; i++)
				{
					_ConvetIndexToSub(i, vuScanLineBase, vuOtherDimLengths);
					size_t uScanLineBase = UConvetSubToIndex(vuScanLineBase, vuDimLengths);
					for(size_t 	b = 0;
								b < UGetNrOfBins(); 
								b++)
							for(size_t 	j = 1, uIndex = uScanLineBase; 
										j < this->vuDimLengths[d]; 
										j++, uIndex += uOffset)
							this->vvBinSATs[b][uIndex + uOffset] += this->vvBinSATs[b][uIndex ];
				}
			}
		}

		//! Allocate the space to store coefficients for all bins. 
		/*! 
		*/
		virtual 
		void 
		_Allocate
		(
			void *_Reserved = NULL
		)
		{
			if( uNrOfBins * uDataSize > uSizeOfFullArrays )
			{
				LOG_ERROR(cerr<<"Exceed the specified application-side memory capacity.");
				exit(EXIT_FAILURE);
			}
			this->vvBinSATs.resize(uNrOfBins);
			for(size_t b = 0; b < uNrOfBins; b++)
				this->vvBinSATs[b].resize(uDataSize);
		}
		
		//! Compute and display statistics for the computed wavelet coefficients.
		virtual 
		void
		_ShowStatistics
		(
			void *_Reserved = NULL
		)
		{
		}

		//! Return the sum of all bins at the given position
		virtual
		void
		_GetAllSums
		(
			const vector<size_t>& vuPos,
			vector<ST>& vdSums,
			void *_Reserved = NULL
		)
		{
			size_t uIndex = UConvetSubToIndex(vuPos, vuDimLengths);
			vdSums.clear();
			for(size_t b = 0; b < UGetNrOfBins(); b++)
				vdSums.push_back(this->vvBinSATs[b][uIndex]);
		}		
		
		CSATEncoder()
		{
		}
	};
}
