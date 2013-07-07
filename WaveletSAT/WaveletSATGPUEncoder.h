#pragma once

#if	WITH_CUDA
#include "CudaDWT.h"
using namespace CudaDWT;
#endif	// #if	WITH_CUDA

#include "libclock.h"	// ADD-BY-LEETEN 2013/07/06
#include "WaveletSATEncoder.h"	// ADD-BY-LEETEN 12/16/2012

namespace WaveletSAT
{
	//! The class that apply WaveletSAT on GPUs
	template<
		typename DT,					//!< Type of the data
		typename ST = double,			//!< Type of the sum
		typename BT = unsigned short,	//!< Type of the bin
		typename WT = double			//!< Type of the wavelet coefficientsd
	>
	class CWaveletSATGPUEncoder:
		#if	WITH_CUDA
		virtual public CCudaDWT,
		#endif
		virtual public CWaveletSATEncoder<DT, ST, BT, WT>	
	{
protected:	
		bool	bIsUsingGPUs;

		int		iTimingPrintingLevel;

		size_t uMaxNrOfElementsOnTheDevice;

		size_t uNrOfElements;

		vector<uint4>			vu4BinSubs;
		vector<typeWavelet>			vfWeights;
		vector<unsigned int>	vuKeys;
		vector<typeWavelet>			vfCoefs;

		vector<unsigned int>	vuSegCounts;	// ADD-BY-LEETEN 01/13/2013

		virtual
		void
		_UpdateBinsOnGPUs
		(
			void *_Reserved = NULL
		)
		{
			// for each wavelet, find the corresponding subscripts
			unsigned int puLevels[CudaDWT::GPU_MAX_NR_OF_DIMS];
			unsigned int puWaveletLengths[CudaDWT::GPU_MAX_NR_OF_DIMS];
				
			// ADD-BY-LEETEN 01/18/2012-BEGIN
			unsigned int puCoefLengths[CudaDWT::GPU_MAX_NR_OF_DIMS];
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				puCoefLengths[d] = (unsigned int)vuCoefLengths[d];
			// ADD-BY-LEETEN 01/18/2012-END
			CCudaDWT::_InitEncoder(
				// ADD-BY-LEETEN 01/18/2012-BEGIN
				UGetNrOfDims(), 
				&puCoefLengths[0],
				// ADD-BY-LEETEN 01/18/2012-END
				uNrOfElements,
				vu4BinSubs.data(),
				vfWeights.data()
				);

			for(size_t c = 0; c < uNrOfUpdatingCoefs; c++)
			{
				bool bIsPrintingTiming = ( iTimingPrintingLevel > 0 )?true:false;	// ADD-BY-LEETEN 01/11/2013
				LIBCLOCK_INIT(bIsPrintingTiming, __FUNCTION__);

				LIBCLOCK_BEGIN(bIsPrintingTiming);
				vector<size_t> vuLevels;
				_ConvertIndexToSub(c, vuLevels, vuDimLevels);
				for(size_t d = 0; d < UGetNrOfDims(); d++)
				{
					puLevels[d] = (unsigned int)vuLevels[d];
					puWaveletLengths[d] = 1 << (unsigned int)(vuDimLevels[d] - vuLevels[d]);
				}
				LIBCLOCK_END(bIsPrintingTiming);

				LIBCLOCK_BEGIN(bIsPrintingTiming);
				size_t uNrOfEncodedCoefs;
				CCudaDWT::_Encode(
					uNrOfElements,
					UGetNrOfDims(),
					&puLevels[0],
					&puWaveletLengths[0],

					&uNrOfEncodedCoefs,
					vuKeys.data(),
					vfCoefs.data(),

					vuSegCounts.data(),	// ADD-BY-LEETEN 01/13/2013

					iTimingPrintingLevel - 1
				);
				LIBCLOCK_END(bIsPrintingTiming);

				LIBCLOCK_BEGIN(bIsPrintingTiming);
				vector<size_t> vuPos;
				vuPos.resize(UGetNrOfDims());
				vector<size_t> vuCoefHalfLengths;
				vuCoefHalfLengths.resize(UGetNrOfDims());
				for(size_t d = 0; d < UGetNrOfDims(); d++)
					vuCoefHalfLengths[d] = vuCoefLengths[d]/2;
				for(size_t e = 0; e < uNrOfEncodedCoefs; e++)
				{
					size_t uKey = (size_t)vuKeys[e];
					unsigned int uCount = vuSegCounts[e];	// ADD-BY-LEETEN 01/13/2013
					for(size_t d = UGetNrOfDims(); d > 0; d--)
					{
						size_t uCoefHalfLength = vuCoefHalfLengths[d - 1];
						vuPos[d - 1] = uKey % uCoefHalfLength;
						uKey /= uCoefHalfLength;
					}
					BT uBin = (BT)uKey;
					this->vcCoefPools[c]._AddAt(uBin, vuPos, (WT)vfCoefs[e], (size_t)uCount);
				}
				LIBCLOCK_END(bIsPrintingTiming);

				LIBCLOCK_PRINT(bIsPrintingTiming);
			}
			uNrOfElements = 0;
		}

		//! Update the specified bin.
		virtual 
		void 
		_UpdateBin
		(
			const vector<size_t>& vuPos, 
			const DT& value,
			const BT& uBin, 
			const WT& dWeight,
			void *_Reserved = NULL
		)
		{
			if( !bIsUsingGPUs )
			{
				CWaveletSATEncoder<DT, ST, BT, WT>::_UpdateBin(vuPos, value, uBin, dWeight);
				return;
			}

			// if the pool for GPU is not full yet, append the tuples of <bin, value, subscripts> to the pool
			vu4BinSubs[uNrOfElements].x = (unsigned int)uBin;
			unsigned int* puBinSub = &vu4BinSubs[uNrOfElements].y;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
				puBinSub[d] = (unsigned int)vuPos[d];
			vfWeights[uNrOfElements] = (typeWavelet)dWeight;
			uNrOfElements++;

			// otherwise, 
			if( uNrOfElements == uMaxNrOfElementsOnTheDevice )
			{
				_UpdateBinsOnGPUs();
			}
		}
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0E00,
			TIMING_PRINTING_LEVEL,	// ADD-BY-LEETEN 01/11/2013
			IS_USING_GPUS,
			MAX_NR_OF_ELEMENTS_ON_THE_DEVICE,
			PARAMETER_END
		};

		virtual	
		void
		_SetInteger(
			int eName,
			long lValue,
			void* _Reserved = NULL
		)
		{
			CWaveletSATEncoder<DT, ST, BT, WT>::_SetInteger(eName, lValue);
			switch(eName)
			{
			// ADD-BY-LEETEN 01/11/2013-BEGIN
			case TIMING_PRINTING_LEVEL:
				iTimingPrintingLevel = lValue;
				break;
			// ADD-BY-LEETEN 01/11/2013-END
			case IS_USING_GPUS:
				bIsUsingGPUs = (lValue)?true:false;
				break;

			case MAX_NR_OF_ELEMENTS_ON_THE_DEVICE:
				uMaxNrOfElementsOnTheDevice = (unsigned int)lValue;
				break;
			}
		}

		//! Finalize the computation of SAT
		virtual	
		void 
		_Finalize
		(
			void *_Reserved = NULL
		)
		{
			if( bIsUsingGPUs && uNrOfElements > 0 )
			{
				// clean up the remaining coefficients in the pool
				_UpdateBinsOnGPUs();
			}

			// call the default _Finalize()
			CWaveletSATEncoder<DT, ST, BT, WT>::_Finalize();
		}

		//! Allocate the space on the GPU device and setup CUDPP
		/*! 
		*/
		virtual	// ADD-BY-LEETEN 09/29/2012
		void 
		_Allocate
		(
			void *_Reserved = NULL
		)
		{
			CWaveletSATEncoder<DT, ST, BT, WT>::_Allocate();
			if( UGetNrOfDims() > CudaDWT::GPU_MAX_NR_OF_DIMS )
				bIsUsingGPUs = false;

			if( !bIsUsingGPUs )
				return;

			uMaxNrOfElementsOnTheDevice = min(uMaxNrOfElementsOnTheDevice, this->uDataSize);
			CCudaDWT::_Init(&uMaxNrOfElementsOnTheDevice);
			vu4BinSubs.resize(uMaxNrOfElementsOnTheDevice);
			vfWeights.resize(uMaxNrOfElementsOnTheDevice);
			vuKeys.resize(uMaxNrOfElementsOnTheDevice);
			vfCoefs.resize(uMaxNrOfElementsOnTheDevice);
			vuSegCounts.resize(uMaxNrOfElementsOnTheDevice);	// ADD-BY-LEETEN 01/13/2013
		}
		
		CWaveletSATGPUEncoder():
			bIsUsingGPUs(false),
			CWaveletSATEncoder<DT, ST, BT, WT>(),
			uNrOfElements(0),
			uMaxNrOfElementsOnTheDevice(CudaDWT::DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE)
		{
		}

		virtual 
			~CWaveletSATGPUEncoder()
		{
			/*
			if( bIsUsingGPUs )
				CudaDWT::_Free();
			*/
		}
	};
}
