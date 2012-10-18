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
	class CBase
	{
		//! #Dimensions
		/*!
		D
		*/
		virtual 
		size_t 
		UGetNrOfDims
		(
			void *_Reserved = NULL	// ADD-BY-LEETEN 10/17/2012
		) const = 0;
		
		//! #Dimensions
		/*!
		B
		*/
		virtual 
		size_t 
		UGetNrOfBins
		(
			void *_Reserved = NULL	// ADD-BY-LEETEN 10/17/2012
		) const = 0;

protected:	
		virtual	
		void 
		_Update
		(
			const vector<size_t>& vuPos,
			const T& value,
			void *_Reserved = NULL
		) = 0;

		//! This method should be overloaded to return the indices of bins that will be changed
		virtual 
		void 
		_MapValueToBins
		(
			const vector<size_t>& vuPos,
			const T& value, 
			vector< pair<size_t, double> >& vuBins,
			void *_Reserved = NULL
		) = 0;
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			NR_OF_ENUMS
		};

		void
		_SetLong(
			int eName,
			long lValue,
			void* _Reserved = NULL
		)
		{
		}

		
		//! Finalize the computation of SAT
		virtual 
		void 
		_Finalize
		(
			void *_Reserved = NULL
		) = 0;

		//! Compute and display statistics for the computed wavelet coefficients.
		virtual 
		void
		_ShowStatistics
		(
		) = 0;

		//! Compute statistics of the compressed result.
		virtual
		void
		_SetDimLengths
		(
			const vector<size_t>& vuDimLengths,
			void *_Reserved = NULL
		) = 0;
		
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
		) = 0;
		
		//! Return the sum of all bins at the given position
		virtual
		void
		_GetAllSums
		(
			const vector<size_t>& vuPos,
			vector<double>& vdSums,
			void *_Reserved = NULL
		) = 0;
	};
}
