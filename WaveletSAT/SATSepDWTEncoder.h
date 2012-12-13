#pragma once

#include <iostream>
#include <vector>
using namespace std;
#include <math.h>

#include "liblog.h"

#include "SATEncoder.h"
#include "SepDWTHeader.h"
#include "WaveletSATEncoder.h"

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

	//! The encoders that computes SAT and applys Sep. DWT in core.
	/*!
	*/
	template<typename DT, typename ST = double>
	class CSATSepDWTEncoder:
		public CSATEncoder<DT, ST>,
		public CWaveletSATEncoder<DT>
	{
protected:	
		////////////////////////////////////////////////////////////////////
		/*
		The protected interface. 
		*/

		virtual	
		void 
		_Update
		(
			const vector<size_t>& vuPos,
			const DT& value,
			void *_Reserved = NULL
		)
		{
			CSATEncoder<DT, ST>::_Update(vuPos, value);
		}

public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/
		enum EParameter
		{
			PARAMETER_BEGIN = 0x0500,
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
//			if( CSATEncoder<DT, ST>::PARAMETER_BEGIN <= eName && eName < CSATEncoder<DT, ST>::PARAMETER_END )
			CSATEncoder<DT, ST>::_SetLong(eName, lValue);

//			if( CWaveletSATEncoder<DT>::PARAMETER_BEGIN <= eName && eName < CWaveletSATEncoder<DT>::PARAMETER_END )
			CWaveletSATEncoder<DT>::_SetLong(eName, lValue);
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
			CSATEncoder<DT, ST>::_Set(vuDimLengths, uNrOfBins);
			CWaveletSATEncoder<DT>::_Set(vuDimLengths, uNrOfBins);
		}

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
			CSATEncoder<DT, ST>::_Finalize();

			CWaveletSATEncoder<DT>::_Allocate();

			size_t uNrOfDims = CWaveletSATEncoder<DT>::UGetNrOfDims();
			size_t uNrOfBins = CWaveletSATEncoder<DT>::UGetNrOfBins();

			vector<size_t> vuSub;
			vuSub.resize(uNrOfDims);
			for(size_t i = 0; i < CSATEncoder<DT, ST>::uDataSize; i++)
			{
				_ConvertIndexToSub(i, vuSub, CSATEncoder<DT, ST>::vuDimLengths);
				size_t uWaveletIndex = UConvertSubToIndex(vuSub, CWaveletSATEncoder<DT>::vuDimLengths);

				for(size_t b = 0; b < uNrOfBins; b++)
					CWaveletSATEncoder<DT>::vvdBinCoefs[b][uWaveletIndex] = CSATEncoder<DT, ST>::vvBinSATs[b][i];
			}

			vector<size_t> vuScanLineBase;
			vuScanLineBase.resize(uNrOfDims);

			vector<size_t> vuOtherDimLengths;
			vuOtherDimLengths.resize(uNrOfDims);

			for(size_t uOffset = 1, d = 0, uDimLength = CWaveletSATEncoder<DT>::vuDimLengths[d]; 
				d < uNrOfDims; 
				uOffset *= uDimLength, d++)
			{
				vector<ST> vSrc;
				vSrc.resize(uDimLength);

				vector<ST> vDst;
				vDst.resize(uDimLength);

				vector<size_t> vuScanLineIndices;
				vuScanLineIndices.resize(uDimLength);

				size_t uNrOfLevels = CWaveletSATEncoder<DT>::vuDimLevels[d] - 1;

				size_t uNrOfScanLines = CWaveletSATEncoder<DT>::uDataSize / uDimLength;

				vuOtherDimLengths = CWaveletSATEncoder<DT>::vuDimLengths;
				vuOtherDimLengths[d] = 1;

				for(size_t i = 0; i < uNrOfScanLines; i++)
				{
					_ConvertIndexToSub(i, vuScanLineBase, vuOtherDimLengths);
					vuScanLineIndices[0] = UConvertSubToIndex(vuScanLineBase, CWaveletSATEncoder<DT>::vuDimLengths);

					// store the element indices along the current scanline so it can be reused for all bins
					for(size_t j = 1; j < uDimLength; j++)
						vuScanLineIndices[j] = vuScanLineIndices[j - 1] + uOffset;

					for(size_t b = 0; b < uNrOfBins; b++)
					{
						for(size_t j = 0; j < uDimLength; j++)
							vSrc[j] = CWaveletSATEncoder<DT>::vvdBinCoefs[b][vuScanLineIndices[j]];

						_DWT1D(vSrc, vDst, uDimLength, uNrOfLevels - 1);

						for(size_t j = 0; j < uDimLength; j++)
							CWaveletSATEncoder<DT>::vvdBinCoefs[b][vuScanLineIndices[j]] = vDst[j];
					}
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
			CSATEncoder<DT, ST>::_Allocate();
		}
		
		//! Compute and display statistics for the computed wavelet coefficients.
		virtual 
		void
		_ShowStatistics
		(
			void *_Reserved = NULL
		)
		{
			CWaveletSATEncoder<DT>::_ShowStatistics();
		}

		//! Retunr the smallest value. 
		/*!
		This can be used to filter too small value caused by numerical error.
		*/
		virtual	
		double 
		DGetThreshold
		(
			void *_Reserved = NULL
		)
		{
			return CWaveletSATEncoder<DT>::DGetThreshold();
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
			CWaveletSATEncoder<DT>::_GetAllSums(vuPos, vdSums);
		}		
		
		CSATSepDWTEncoder():
			CSATEncoder<DT, ST>(),
			CWaveletSATEncoder<DT>()
		{
		}
	};
}
