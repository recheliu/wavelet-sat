#pragma once

#if	WITH_CUDA
#include "CudaDWT.h"
using namespace CudaDWT;
#endif	// #if	WITH_CUDA

#define	WITH_CPU_BUCKET_SORT	0	// ADD-BY-LEETEN 01/13/2013

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

		// MOD-BY-LEETEN 01/11/2013-FROM:		bool	bIsPrintingTiming;
		int		iTimingPrintingLevel;
		// MOD-BY-LEETEN 01/11/2013-END

		size_t uMaxNrOfElementsOnTheDevice;

		size_t uNrOfElements;

		vector<uint4>			vu4BinSubs;
		vector<float>			vfWeights;
		vector<unsigned int>	vuKeys;
		vector<float>			vfCoefs;

		vector<unsigned int>	vuSegCounts;	// ADD-BY-LEETEN 01/13/2013

		// ADD-BY-LEETEN 01/13/2013-BEGIN
		#if		WITH_CPU_BUCKET_SORT
		vector<uint4>			vu4SortedBinSubs;
		vector<float>			vfSortedWeights;

		vector<size_t>			vuBucketCounts;

		vector<size_t>			vuBucketBases;
		#endif	// #if	WITH_CPU_BUCKET_SORT
		// ADD-BY-LEETEN 01/13/2013-END

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
				
			#if		!WITH_CPU_BUCKET_SORT	// ADD-BY-LEETEN 01/13/2013
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
			// ADD-BY-LEETEN 01/13/2013-BEGIN
			#else	// #if	!WITH_CPU_BUCKET_SORT	
			vuBucketBases[0] = 0;
			for(size_t b = 1; b < vuBucketCounts.size(); b++)
				vuBucketBases[b] = vuBucketBases[b - 1] + vuBucketCounts[b - 1];
			fill(vuBucketCounts.begin(), vuBucketCounts.end(), 0);	// reuse this buffer to hold the count
			for(size_t d = 0; d < vu4BinSubs.size(); d++)
			{
				size_t uBin = vu4BinSubs[d].x;
				size_t uSortedIndex = vuBucketBases[uBin] + vuBucketCounts[uBin];
				vu4SortedBinSubs[uSortedIndex] = vu4BinSubs[d];
				vfSortedWeights[uSortedIndex] = vfWeights[d];
				vuBucketCounts[uBin]++;
			}
			// now reset the counts to 0
			fill(vuBucketCounts.begin(), vuBucketCounts.end(), 0);

			CCudaDWT::_InitEncoder(
				uNrOfElements,
				vu4SortedBinSubs.data(),
				vfSortedWeights.data(),
				true
				);
			#endif	// #if	!WITH_CPU_BUCKET_SORT	
			// ADD-BY-LEETEN 01/13/2013-END

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
					// MOD-BY-LEETEN 01/11/2013-FROM:	vfCoefs.data()
					vfCoefs.data(),

					vuSegCounts.data(),	// ADD-BY-LEETEN 01/13/2013

					iTimingPrintingLevel - 1
					// MOD-BY-LEETEN 01/11/2013-END
				);
				LIBCLOCK_END(bIsPrintingTiming);

				LIBCLOCK_BEGIN(bIsPrintingTiming);
				for(size_t e = 0; e < uNrOfEncodedCoefs; e++)
				{
					vector<size_t> vuPos;
					vuPos.resize(UGetNrOfDims());
					size_t uKey = (size_t)vuKeys[e];
					unsigned int uCount = vuSegCounts[e];	// ADD-BY-LEETEN 01/13/2013
					#if	0	// MOD-BY-LEETEN 01/18/2012-FROM:
					for(size_t d = 0; d < UGetNrOfDims(); d++, uKey /= 256)
						vuPos[UGetNrOfDims() - 1 - d] = uKey % 256;
					#else	// MOD-BY-LEETEN 01/18/2012-TO:
					for(size_t d = UGetNrOfDims(); d > 0; d--)
					{
						size_t uCoefHalfLength = vuCoefLengths[d - 1]/2;
						vuPos[d - 1] = uKey % uCoefHalfLength;
						uKey /= uCoefHalfLength;
					}
					#endif	// MOD-BY-LEETEN 01/18/2012-END
					BT uBin = (BT)uKey;
					// MOD-BY-LEETEN 01/13/2013-FROM:					this->vcCoefPools[c]._AddAt(uBin, vuPos, (WT)vfCoefs[e]);
					this->vcCoefPools[c]._AddAt(uBin, vuPos, (WT)vfCoefs[e], (size_t)uCount);
					// MOD-BY-LEETEN 01/13/2013-END
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

			// ADD-BY-LEETEN 01/13/2013-BEGIN
			#if	WITH_CPU_BUCKET_SORT
			// update the bucket size of this bin
			vuBucketCounts[uBin]++;
			#endif	// #if	WITH_CPU_BUCKET_SORT
			// ADD-BY-LEETEN 01/13/2013-END

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

			// ADD-BY-LEETEN 01/13/2013-BEGIN
			#if		WITH_CPU_BUCKET_SORT
			// allocate the buffer to hold the sorted and consecutive keys
			vu4SortedBinSubs.resize(uMaxNrOfElementsOnTheDevice);
			// allocate the buffer to hold the corresponding values
			vfSortedWeights.resize(uMaxNrOfElementsOnTheDevice);

			// allocate the buffers for bucket sort
			vuBucketBases.assign(UGetNrOfBins(), 0);
			vuBucketCounts.assign(UGetNrOfBins(), 0);
			#endif	// #if	WITH_CPU_BUCKET_SORT
			// ADD-BY-LEETEN 01/13/2013-END
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
