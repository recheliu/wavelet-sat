#include <vector>
using namespace std;
#include <stdio.h>
#include <assert.h>
#include "libclock.h"
#include "liblog.h"
#include "cuda_macro.h"
#include "CudaDWT.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

#define BLOCK_SIZE_X	16
#define BLOCK_SIZE_Y	8
#define BLOCK_SIZE		(BLOCK_SIZE_X * BLOCK_SIZE_Y)

#include "ProjToWavelet_kernel.h"
#include "MarkSegments_kernel.h"

namespace CudaDWT
{
	CCudaDWT::
		~CCudaDWT
	(
	)
	{
		if( bIsInitialized )
		{
			// free the memory space
			FREE_MEMORY(pu4BinSub_device);
			FREE_MEMORY(pfValues_device);
			FREE_MEMORY(puKeys_device);
			FREE_MEMORY(pfCoefs_device);
			FREE_MEMORY(pfCompactedCoefs_device);
			FREE_MEMORY(puCompactedKeys_device);

			FREE_MEMORY(puOnes_device);
			FREE_MEMORY(puCompactedSegCounts_device);
		}
	}

	void
	CCudaDWT::
	_Init
	(
		size_t* puMaxNrOfElementsOnTheDevice,
		void* _Reserved
	)
	{
		if( *puMaxNrOfElementsOnTheDevice > CudaDWT::DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE )
		{
			LOG_ERROR(cerr<<"uMaxNrOfElementsOnTheDevice is clampped to CudaDWT::DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE");
			*puMaxNrOfElementsOnTheDevice = CudaDWT::DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE;
		}

		size_t uMaxNrOfElementsOnTheDevice = *puMaxNrOfElementsOnTheDevice;
		// allocate the memory space
		CUDA_SAFE_CALL(cudaMalloc((void**)&pu4BinSub_device,		sizeof(pu4BinSub_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&pfValues_device,			sizeof(pfValues_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&puKeys_device,			sizeof(puKeys_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&pfCoefs_device,			sizeof(pfCoefs_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&pfCompactedCoefs_device,	sizeof(pfCompactedCoefs_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&puCompactedKeys_device,	sizeof(puCompactedKeys_device[0]) * uMaxNrOfElementsOnTheDevice));

		CUDA_SAFE_CALL(cudaMalloc((void**)&puOnes_device,			sizeof(puOnes_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&puCompactedSegCounts_device,			sizeof(puCompactedSegCounts_device[0]) * uMaxNrOfElementsOnTheDevice));
		vector<unsigned int> vuOnes;
		vuOnes.assign(uMaxNrOfElementsOnTheDevice, 1);
		CUDA_SAFE_CALL(cudaMemcpy(puOnes_device, vuOnes.data(), sizeof(puOnes_device[0]) * vuOnes.size(), cudaMemcpyHostToDevice));

		bIsInitialized = true;
	}

	void
	CCudaDWT::
	_InitEncoder
	(
		size_t uNrOfDims,
		unsigned int puCoefLengths[],
		size_t				uNrOfElements,
		const uint4			pu4BinSubs[],
		const typeValue	pfValues[],
		bool bWithCpuBucketSort,	
		void* _Reserved
	)
	{
		this->bWithCpuBucketSort = bWithCpuBucketSort;	

		// upload the tuples in the pool to the device side
		CUDA_SAFE_CALL(
			cudaMemcpy(
				&pu4BinSub_device[0], 
				&pu4BinSubs[0],		
				sizeof(pu4BinSub_device[0]) * uNrOfElements, 
				cudaMemcpyHostToDevice));

		CUDA_SAFE_CALL(
			cudaMemcpy(
				&pfValues_device[0],	
				&pfValues[0],
				sizeof(pfValues_device[0]) * uNrOfElements, 
				cudaMemcpyHostToDevice));

		CUDA_SAFE_CALL(
			cudaMemcpyToSymbol(
				puCoefLengths_const, 
				puCoefLengths,
				uNrOfDims * sizeof(puCoefLengths_const[0]), 
				0,
				cudaMemcpyHostToDevice));

		v3Blk = dim3(BLOCK_SIZE);
		size_t uNrOfBlocks = (size_t)ceilf((float)uNrOfElements / (float)v3Blk.x);
		size_t uGridSizeX = (size_t)ceil(sqrtf((float)uNrOfBlocks));
		size_t uGridSizeY = (size_t)ceil((float)uNrOfBlocks/(float)uGridSizeX);
		v3Grid = dim3( (unsigned int)uGridSizeX, (unsigned int)uGridSizeY );
	}

	void
	CCudaDWT::
	_Encode
	(
		size_t				uNrOfBins,	
		size_t				uNrOfElements,
		size_t				uNrOfDims,
		const unsigned int	puLevels[],
		const unsigned int	puWaveletLengths[],

		size_t				*puNrOfElements,
		typeKey				puKeys_host[],
		typeCoef			pfCoefs_host[],
		unsigned int		puSegCounts_host[],

		int iTimingPrintingLevel,	
		void* _Reserved
	)
	{
		bool bIsPrintingTiming = (iTimingPrintingLevel > 0)?true:false;	
		LIBCLOCK_INIT(bIsPrintingTiming, __FUNCTION__);

		// copy the lengths of the local coefficient array, wavelet lengths, and levels
		LIBCLOCK_BEGIN(bIsPrintingTiming);

		CUDA_SAFE_CALL(
			cudaMemcpyToSymbol(
				puLevels_const, 
				&puLevels[0], 
				sizeof(puLevels_const[0]) * uNrOfDims,
				0,
				cudaMemcpyHostToDevice) );

		CUDA_SAFE_CALL(
			cudaMemcpyToSymbol(
				puWaveletLengths_const, 
				&puWaveletLengths[0], 
				sizeof(puWaveletLengths_const[0]) * uNrOfDims, 
				0,
				cudaMemcpyHostToDevice) );
		LIBCLOCK_END(bIsPrintingTiming);

		LIBCLOCK_BEGIN(bIsPrintingTiming);
		// 

		_ProjToWavelet_kernel<<<v3Grid, v3Blk, 0>>>(
			&pu4BinSub_device[0],	// the tuples of <bin, data_subscripts> of all elements
			&pfValues_device[0],	// the counts of all elements
			(unsigned int)uNrOfBins,		
			(unsigned int)uNrOfDims, 
			(unsigned int)uNrOfElements,
			&puKeys_device[0],		// output: the keys of all elements. The keys are composed of bin and local_subscripts
			&pfCoefs_device[0]		// output: the wavelet projection of the current wavelet 
			);
		CUT_CHECK_ERROR("_ProjToWavelet_kernel() failed");
		LIBCLOCK_END(bIsPrintingTiming);

		if( !bWithCpuBucketSort )	
		{
			// sort the wavelet projection according to the key composed by the bin and local subscripts
			LIBCLOCK_BEGIN(bIsPrintingTiming);
			thrust::device_ptr<typeKey> vKeys_device(puKeys_device);
			thrust::device_ptr<typeCoef> vCoefs_device(pfCoefs_device);
			thrust::sort_by_key(vKeys_device, vKeys_device + uNrOfElements, vCoefs_device);
			LIBCLOCK_END(bIsPrintingTiming);
		}	

		thrust::device_ptr<typeKey> vuKeys_device(puKeys_device);
		thrust::device_ptr<unsigned int> vuOnes_device(puOnes_device);
		thrust::device_ptr<typeKey> vuCompactedKeys_device(puCompactedKeys_device);
		thrust::device_ptr<unsigned int> vuCompactedSegCounts_device(puCompactedSegCounts_device);
		thrust::pair< thrust::device_ptr<typeKey>, thrust::device_ptr<unsigned int> > pairEnds = 
			thrust::reduce_by_key<
				thrust::device_ptr<typeKey>, 
				thrust::device_ptr<unsigned int>, 
				thrust::device_ptr<typeKey>, 
				thrust::device_ptr<unsigned int> >
			(
				vuKeys_device, 
				vuKeys_device + uNrOfElements, 
				vuOnes_device,
				vuCompactedKeys_device,
				vuCompactedSegCounts_device
				);
		size_t uNrOfCompactedSegCounts_host = size_t(pairEnds.second - vuCompactedSegCounts_device);
		size_t uNrOfCompactedKeys_host = size_t(pairEnds.first - vuCompactedKeys_device);

		LIBCLOCK_BEGIN(bIsPrintingTiming);
		CUDA_SAFE_CALL(
			cudaMemcpy(
				&puSegCounts_host[0],
				&puCompactedSegCounts_device[0], 
				uNrOfCompactedSegCounts_host * sizeof(puSegCounts_host[0]),
				cudaMemcpyDeviceToHost) );
		LIBCLOCK_END(bIsPrintingTiming);
		thrust::device_ptr<typeCoef> vfCoefs_device(pfCoefs_device);
		thrust::device_ptr<typeCoef> vfCompactedCoefs_device(pfCompactedCoefs_device);
		size_t uNrOfCompactedCoefs_host;
		{
			thrust::pair<thrust::device_ptr<typeKey>, thrust::device_ptr<typeCoef> > pairEnds = 
				thrust::reduce_by_key<
					thrust::device_ptr<typeKey>, 
					thrust::device_ptr<typeCoef>, 
					thrust::device_ptr<typeKey>, 
					thrust::device_ptr<typeCoef> >
				(
					vuKeys_device, 
					vuKeys_device + uNrOfElements, 
					vfCoefs_device,
					vuCompactedKeys_device,
					vfCompactedCoefs_device);
			uNrOfCompactedCoefs_host = size_t(pairEnds.second - vfCompactedCoefs_device);
		}
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		CUDA_SAFE_CALL(
			cudaMemcpy(
				&pfCoefs_host[0],
				&pfCompactedCoefs_device[0], 
				uNrOfCompactedCoefs_host * sizeof(pfCoefs_host[0]),
				cudaMemcpyDeviceToHost) );
		LIBCLOCK_END(bIsPrintingTiming);

		// download the keys and the coefficinets back 
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		CUDA_SAFE_CALL(
			cudaMemcpy(
				&puKeys_host[0],
				&puCompactedKeys_device[0], 
				uNrOfCompactedKeys_host * sizeof(puKeys_host[0]),
				cudaMemcpyDeviceToHost) );
		LIBCLOCK_END(bIsPrintingTiming);

		ASSERT_OR_LOG(uNrOfCompactedKeys_host == uNrOfCompactedCoefs_host, cerr<<"Unmatched #keys and #coefs.");
		*puNrOfElements = uNrOfCompactedKeys_host;
		LIBCLOCK_PRINT(bIsPrintingTiming);
	}
};


/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
