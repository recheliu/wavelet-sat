#pragma once

#include <vector>
using namespace std;
#include <math.h>

#include "HeaderBase.h"

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

	//! The base class of the encoders
	/*!
	In order to apply wavelet transform for SATs or Integral Histograms, please follow the procedure below:
	1. Setup the data size by _SetDimLengths().
	2. Setup the #SATs by _AllocateBins().
	3. Call _Update() for each data point to incrementally compute wavelet.
	4. Call _Finalize() to finish the computation.

	To query the entry for a certain location:
	1. Call _GetAllSums() to get all SAT values for the given location.
	*/
	// MOD-BY-LEETEN 01/03/2013-FROM:	template<typename DT, typename ST>
	template<
		typename DT,	//!< Type of the data
		typename ST = typeSum,	//!< Type of the sum
		typename BT = typeBin	//!< Type of the bin
	>
	// MOD-BY-LEETEN 01/03/2013-END
	class CEncoderBase
		:virtual public CBase	// ADD-BY-LEETEN 12/30/2012
	{
protected:	
		//! Update the specified bin.
		virtual
		void 
		_UpdateBin
		(
			const vector<size_t>& vuPos, 
			const DT& value,
			// MOD-BY-LEETEN 01/03/2013-FROM:			size_t uBin, 
			const BT& uBin, 
			// MOD-BY-LEETEN 01/03/2013-END
			const ST& weight,
			void *_Reserved = NULL
		) = 0;

		virtual 
		void 
		_MapValueToBins
		(
			const vector<size_t>& vuPos,
			const DT& value, 
			// MOD-BY-LEETEN 01/03/2013-FROM:			vector< pair<size_t, ST> >& vpBins,
			vector< pair<BT, ST> >& vpBins,
			// MOD-BY-LEETEN 01/03/2013-END
			void *_Reserved = NULL
		) = 0;

		virtual	
		void 
		_Update
		(
			const vector<size_t>& vuPos,
			const DT& value,
			void *_Reserved = NULL
		)
		{
			// MOD-BY-LEETEN 01/03/2013-FROM:			vector< pair<size_t, ST> > vpBins;
			vector< pair<BT, ST> > vpBins;
			// MOD-BY-LEETEN 01/03/2013-END
			_MapValueToBins(vuPos, value, vpBins);
			// MOD-BY-LEETEN 01/03/2013-FROM:			for(typename vector< pair<size_t, ST> >::iterator ivpBin  = vpBins.begin(); ivpBin != vpBins.end(); ivpBin++)
			for(typename vector< pair<BT, ST> >::iterator ivpBin  = vpBins.begin(); ivpBin != vpBins.end(); ivpBin++)
			// MOD-BY-LEETEN 01/03/2013-END
				_UpdateBin(vuPos, value, ivpBin->first, ivpBin->second);
		}

	public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/
		enum EParameter
		{
			PARAMETER_BEGIN = 0x0200,
			PARAMETER_END
		};

		//! Finalize the computation of SAT
		virtual 
		void 
		_Finalize
		(
			void *_Reserved = NULL
		) = 0;

		//! Allocate the space to store/compute SATs
		virtual 
		void 
		_Allocate
		(
			void *_Reserved = NULL
		) = 0;

		//! Compute and display statistics for the computed wavelet coefficients.
		virtual 
		void
		_ShowStatistics
		(
			void *_Reserved = NULL
		) = 0;

		// ADD-BY-LEETEN 12/12/2012-BEGIN
		//! Save the coefficients to a file
		virtual 
		void
		_SaveFile
		(
		 const char* szFilepathPrefix, 
			void *_Reserved = NULL
		) = 0;
		// ADD-BY-LEETEN 12/12/2012-END
	};
}
