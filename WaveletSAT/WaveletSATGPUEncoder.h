#pragma once

#if	WITH_CUDA
#include "CudaDWT.h"
using namespace CudaDWT;
#endif	// #if	WITH_CUDA

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

		bool	bIsPrintingTiming;

		size_t uMaxNrOfElementsOnTheDevice;

		size_t uNrOfElements;

		vector<uint4>			vu4BinSubs;
		vector<float>			vfWeights;
		vector<unsigned int>	vuKeys;
		vector<float>			vfCoefs;

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
				
			CCudaDWT::_InitEncoder(
				uNrOfElements,
				vu4BinSubs.data(),
				vfWeights.data()
				);

			for(size_t c = 0; c < uNrOfUpdatingCoefs; c++)
			{
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
					vfCoefs.data()
				);
				LIBCLOCK_END(bIsPrintingTiming);

				LIBCLOCK_BEGIN(bIsPrintingTiming);
				for(size_t e = 0; e < uNrOfEncodedCoefs; e++)
				{
					vector<size_t> vuPos;
					vuPos.resize(UGetNrOfDims());
					size_t uKey = (size_t)vuKeys[e];
					for(size_t d = 0; d < UGetNrOfDims(); d++, uKey /= 256)
						vuPos[UGetNrOfDims() - 1 - d] = uKey % 256;
					BT uBin = (BT)uKey;
					this->vcCoefPools[c]._AddAt(uBin, vuPos, (WT)vfCoefs[e]);
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
			vfWeights[uNrOfElements] = (float)dWeight;
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
