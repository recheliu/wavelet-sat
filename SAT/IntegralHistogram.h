#pragma once

#include <vector>
using namespace std;
#include <math.h>

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace SAT
{
	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

	//! The base class of WaveletSAT.
	/*!
	In order to apply wavelet transform for SATs or Integral Histograms, please follow the procedure below:
	1. Setup the data size by _SetDimLengths().
	2. Setup the #SATs by _AllocateBins().
	3. Call _Update() for each data point to incrementally compute wavelet.
	4. Call _Finalize() to finish the computation.

	To query the entry for a certain location:
	1. Call _GetAllSums() to get all SAT values for the given location.
	*/
	template<class T>
	class CIntegralHistogram:
		public CBase<T>
	{
protected:	
		size_t uSizeOfFullArrays;	// ADD-BY-LEETEN 10/18/2012

		size_t uDataLength;

		vector<size_t> vuDimLengths;

		vector< vector<T> > vvdBinSATs;

		//! #Dimensions
		/*!
		D
		*/
		virtual 
		size_t 
		UGetNrOfDims
		(
			void *_Reserved = NULL
		) const
		{
			return this->vuDimLengths.size();
		}
		
		//! #Dimensions
		/*!
		B
		*/
		virtual 
		size_t 
		UGetNrOfBins
		(
			void *_Reserved = NULL
		) const
		{
			return this->vvdBinSATs.size();
		}

		void
		_ConvetIndexToSub
		(
			size_t uIndex,
			vector<size_t>& vuSub,
			void* _Reserved = NULL
		)
		{
			if( this->UGetNrOfDims() != vuSub.size() )
				vuSub.resize(this->UGetNrOfDims());
			for(size_t d = 0; d < this->UGetNrOfDims(); d++)
			{
				vuSub[d] = uIndex % this->vuDimLengths[d];
				uIndex /= this->vuDimLengths[d];
			}
		};

		size_t
		UConvetSubToIndex
		(
			const vector<size_t>& vuSub,
			void* _Reserved = NULL
		)
		{
			size_t uIndex = 0;
			for(size_t d = 0, vuSubSize = 1; d < vuSub.size(); vuSubSize *= this->vuDimLengths[d], d++)
				uIndex += vuSubSize * vuSub[d];
			return uIndex;
		}

		//! Update the specified bin.
		void 
		_UpdateBin
		(
			const vector<size_t>& vuPos, 
			const T& value,
			size_t uBin, 
			double dWeight,
			void *_Reserved = NULL
		)
		{
			size_t uIndex = this->UConvetSubToIndex(vuPos);
			this->vvdBinSATs[uBin][uIndex] += dWeight;
		}

		virtual	
		void 
		_Update
		(
			const vector<size_t>& vuPos,
			const T& value,
			void *_Reserved = NULL
		)
		{
			vector< pair<size_t, double> > vpBins;

			_MapValueToBins(vuPos, value, vpBins);

			for(vector< pair<size_t, double> >::iterator ivpBin  = vpBins.begin(); ivpBin != vpBins.end(); ivpBin++)
				_UpdateBin(vuPos, value, ivpBin->first, ivpBin->second);
		}

public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			//! Total size in MB of all bin SATs.
			SIZE_OF_FULL_ARRAYS,

			NR_OF_ENUMS
		};

		void
		_SetLong(
			int eName,
			long lValue,
			void* _Reserved = NULL
		)
		{
			// ADD-BY-LEETEN 10/18/2012-BEGIN
			switch(eName)
			{
			case SIZE_OF_FULL_ARRAYS:
				uSizeOfFullArrays = (size_t)lValue * 1024 * 1024;
				break;
			}
			// ADD-BY-LEETEN 10/18/2012-END
		}

		double 
		DGetThreshold
		(
			void *_Reserved = NULL
		)
		{
			return 0.0;
		}
		
		//! Finalize the computation of SAT
		virtual 
		void 
		_Finalize
		(
			void *_Reserved = NULL
		)
		{
			vector<size_t> vuScanLineBase;	// the base coordinate of the scane line
			vuScanLineBase.resize(this->UGetNrOfDims());

			for(size_t d = 0; d < this->UGetNrOfDims(); d++)
			{
				size_t uNrOfScanLines = this->uDataLength / this->vuDimLengths[d];

				vector<size_t> vuScanLineIndices;
				vuScanLineIndices.resize(this->vuDimLengths[d]);

				for(size_t i = 0; i < uNrOfScanLines; i++)
				{
					for(size_t uIndex = i, d2 = 0; d2 < this->UGetNrOfDims(); d2++)
					{
						if( d2 == d )
							vuScanLineBase[d2] = 0;
						else
						{
							vuScanLineBase[d2] = (uIndex % this->vuDimLengths[d2]);
							uIndex /= this->vuDimLengths[d2];
						}
					}

					// store the element indices along the current scanline so it can be reused for all bins
					for(size_t j = 0; j < this->vuDimLengths[d]; j++)
					{
						vuScanLineBase[d] = j;
						vuScanLineIndices[j] = this->UConvetSubToIndex(vuScanLineBase);
					}

					for(size_t b = 0; b < UGetNrOfBins(); b++)
						for(size_t j = 1; j < this->vuDimLengths[d]; j++)
							this->vvdBinSATs[b][vuScanLineIndices[j]] += this->vvdBinSATs[b][vuScanLineIndices[j - 1]];
				}
			}
		}

		//! Compute and display statistics for the computed wavelet coefficients.
		virtual 
		void
		_ShowStatistics
		(
		)
		{
		}

		//! Compute statistics of the compressed result.
		virtual
		void
		_SetDimLengths
		(
			const vector<size_t>& vuDimLengths,
			void *_Reserved = NULL
		)
		{
			this->vuDimLengths.resize(vuDimLengths.size());
			uDataLength = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				this->vuDimLengths[d] = vuDimLengths[d];
				uDataLength *= vuDimLengths[d];
			}
		}
		
		//! Allocate the space to store coefficients for all bins. 
		/*! 
		*/
		virtual 
		void 
		_AllocateBins
		(
			size_t uNrOfBins,
			vector<size_t>& vuDimMaxLevels,
			void *_Reserved = NULL
		)
		{
			// ADD-BY-LEETEN 10/18/2012-BEGIN
			if( uNrOfBins * uDataLength > uSizeOfFullArrays )
			{
				LOG_ERROR(cerr<<"Exceed the specified application-side memory capacity.");
				exit(EXIT_FAILURE);
			}
			// ADD-BY-LEETEN 10/18/2012-END
			this->vvdBinSATs.resize(uNrOfBins);
			for(size_t b = 0; b < uNrOfBins; b++)
				this->vvdBinSATs[b].resize(uDataLength);
		}
		
		//! Return the sum of all bins at the given position
		virtual
		void
		_GetAllSums
		(
			const vector<size_t>& vuPos,
			vector<double>& vdSums,
			void *_Reserved = NULL
		)
		{
			size_t uIndex = this->UConvetSubToIndex(vuPos);
			vdSums.clear();
			for(size_t b = 0; b < UGetNrOfBins(); b++)
				vdSums.push_back(this->vvdBinSATs[b][uIndex]);
		}

		CIntegralHistogram()
		{
		}
	};
}
