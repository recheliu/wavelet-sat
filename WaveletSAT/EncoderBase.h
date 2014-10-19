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
	template<
		typename DT,	//!< Type of the data
		typename ST = typeSum,	//!< Type of the sum
		typename BT = typeBin	//!< Type of the bin
	>
	class CEncoderBase
		:virtual public CBase	
	{
protected:	
		//! Update the specified bin.
		virtual
		void 
		_UpdateBin
		(
			const vector<size_t>& vuPos, 
			const DT& value,
			const BT& uBin, 
			const ST& weight,
			void *_Reserved = NULL
		) = 0;

		virtual 
		void 
		_MapValueToBins
		(
			const vector<size_t>& vuPos,
			const DT& value, 
			vector< pair<BT, ST> >& vpBins,
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
			vector< pair<BT, ST> > vpBins;
			_MapValueToBins(vuPos, value, vpBins);
			for(typename vector< pair<BT, ST> >::iterator ivpBin  = vpBins.begin(); ivpBin != vpBins.end(); ivpBin++)
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

		//! Save the coefficients to a file
		virtual 
		void
		_SaveFile
		(
		 const char* szFilepathPrefix, 
			void *_Reserved = NULL
		) = 0;
	};
}

/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
