#pragma once

#include <vector>
using namespace std;
#include <math.h>

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	void
	_ConvetIndexToSub
	(
		size_t uIndex,
		vector<size_t>& vuSub,
		const vector<size_t>& vuDimLengths,
		void* _Reserved = NULL
	)
	{
		size_t uNrOfDims = vuDimLengths.size();
		if( uNrOfDims != vuSub.size() )
			vuSub.resize(uNrOfDims);
		for(size_t d = 0; d < uNrOfDims; d++)
		{
			vuSub[d] = uIndex % vuDimLengths[d];
			uIndex /= vuDimLengths[d];
		}
	};

	size_t
	UConvetSubToIndex
	(
		const vector<size_t>& vuSub,
		const vector<size_t>& vuDimLengths,
		void* _Reserved = NULL
	)
	{
		size_t uNrOfDims = vuDimLengths.size();
		size_t uIndex = 0;
		for(size_t d = 0, vuSubSize = 1; d < vuSub.size(); vuSubSize *= vuDimLengths[d], d++)
			uIndex += vuSubSize * vuSub[d];
		return uIndex;
	}

	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

	//! The base class of the header information.
	/*!
	In order to apply wavelet transform for SATs or Integral Histograms, please follow the procedure below:
	1. Setup the data size by _SetDimLengths().
	2. Setup the #SATs by _AllocateBins().
	3. Call _Update() for each data point to incrementally compute wavelet.
	4. Call _Finalize() to finish the computation.

	To query the entry for a certain location:
	1. Call _GetAllSums() to get all SAT values for the given location.
	*/
	class CHeaderBase
	{
protected:	
		//! The size of the allocated memory for the coefficients
		size_t uSizeOfFullArrays;	

		//! The size of the coefficients
		size_t uDataSize;

		//! The length of the coefficients per dim.
		vector<size_t> vuDimLengths;

		//! The #bins of this SAT
		size_t uNrOfBins;
		
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
			return uNrOfBins;
		}
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0000,

			//! Total size in MB of all bin SATs.
			SIZE_OF_FULL_ARRAYS,

			PARAMETER_END
		};

		virtual
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

		virtual
		double 
		DGetThreshold
		(
			void *_Reserved = NULL
		)
		{
			return 0.0;
		}
		
		//! Compute statistics of the compressed result.
		virtual
		void
		_Set
		(
			const vector<size_t>& vuDimLengths,
			const size_t uNrOfBins,
			void *_Reserved = NULL
		)
		{
			this->vuDimLengths.resize(vuDimLengths.size());
			uDataSize = 1;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				this->vuDimLengths[d] = vuDimLengths[d];
				uDataSize *= vuDimLengths[d];
			}
			this->uNrOfBins = uNrOfBins;
		}
		
		CHeaderBase()
		{
		}
	};
}
