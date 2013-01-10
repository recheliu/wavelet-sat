#include <vector>
using namespace std;
#include <stdio.h>
#include <assert.h>
#include "libclock.h"
#include "liblog.h"
#include "cuda_macro.h"
#include "CudaDWT.h"

inline const char* SZGetCudppError(const CUDPPResult result)
{
	switch(result)
	{
	case CUDPP_ERROR_INVALID_HANDLE:
		return	"Specified handle (for example, to a plan) is invalid.";
	case CUDPP_ERROR_ILLEGAL_CONFIGURATION:
		return	"Specified configuration is illegal. For example, an invalid or illogical combination of options. ";
	case CUDPP_ERROR_INVALID_PLAN:
		return	"The plan is not configured properly. For example, passing a plan for scan to cudppSegmentedScan. ";
	case CUDPP_ERROR_INSUFFICIENT_RESOURCES:
		return	"The function could not complete due to insufficient resources (typically CUDA device resources such as shared memory) for the specified problem size. ";
	case CUDPP_ERROR_UNKNOWN:
		return	"Unknown or untraceable error.";
	}
	return "No Error";
}

#define ASSERT_CUDPP(call)	\
	{	\
		CUDPPResult result = call;	\
		if( CUDPP_SUCCESS != result )	\
		{	\
			CUT_CHECK_ERROR(# call);	\
			LOG_ERROR(cerr<<SZGetCudppError(result));	\
		}	\
	}	\
	

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
			FREE_MEMORY(puiSegFlags_device);
			FREE_MEMORY(pfCoefs_device);
			FREE_MEMORY(pfCoefSums_device);
			FREE_MEMORY(pfCompactedCoefs_device);
			FREE_MEMORY(puCompactedKeys_device);

			FREE_MEMORY(puNrOfCompactedCoefs_device);
			FREE_MEMORY(puNrOfCompactedKeys_device);

			// free cudpp resources
			if( planSort )
				ASSERT_CUDPP(cudppDestroyPlan(planSort));  
			if( planSegScanCoefs )
				ASSERT_CUDPP(cudppDestroyPlan(planSegScanCoefs));
			if( planCompactCoefs )
				ASSERT_CUDPP(cudppDestroyPlan(planCompactCoefs));
			if( planCompactKeys )
				ASSERT_CUDPP(cudppDestroyPlan(planCompactKeys));
			if( theCudpp )
				ASSERT_CUDPP(cudppDestroy(theCudpp));
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
		ASSERT_CUDPP(cudppCreate(&theCudpp));
		
		if( *puMaxNrOfElementsOnTheDevice > CudaDWT::DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE )
		{
			LOG_ERROR(cerr<<"uMaxNrOfElementsOnTheDevice is clampped to CudaDWT::DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE");
			*puMaxNrOfElementsOnTheDevice = CudaDWT::DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE;
		}

		size_t uMaxNrOfElementsOnTheDevice = *puMaxNrOfElementsOnTheDevice;

		configSort.op = CUDPP_ADD;
		configSort.datatype = CUDPP_UINT;
		configSort.algorithm = CUDPP_SORT_RADIX;
		configSort.options = CUDPP_OPTION_KEY_VALUE_PAIRS;
		ASSERT_CUDPP(cudppPlan(theCudpp, &planSort, configSort, uMaxNrOfElementsOnTheDevice, 1, 0));  

		configSegScanCoefs.op = CUDPP_ADD;
		configSegScanCoefs.datatype = CUDPP_FLOAT;
		configSegScanCoefs.algorithm = CUDPP_SEGMENTED_SCAN;
		configSegScanCoefs.options = CUDPP_OPTION_BACKWARD | CUDPP_OPTION_INCLUSIVE;
		ASSERT_CUDPP(cudppPlan(theCudpp, &planSegScanCoefs, configSegScanCoefs, uMaxNrOfElementsOnTheDevice, 1, 0));

		configCompactCoefs.datatype = CUDPP_FLOAT;
		configCompactCoefs.algorithm = CUDPP_COMPACT;
		configCompactCoefs.options = CUDPP_OPTION_FORWARD;
		ASSERT_CUDPP(cudppPlan(theCudpp, &planCompactCoefs, configCompactCoefs, uMaxNrOfElementsOnTheDevice, 1, 0));

		configCompactKeys.datatype = CUDPP_UINT;
		configCompactKeys.algorithm = CUDPP_COMPACT;
		configCompactKeys.options = CUDPP_OPTION_FORWARD;
		ASSERT_CUDPP(cudppPlan(theCudpp, &planCompactKeys, configCompactKeys, uMaxNrOfElementsOnTheDevice, 1, 0));

		// allocate the memory space
		CUDA_SAFE_CALL(cudaMalloc((void**)&pu4BinSub_device,		sizeof(pu4BinSub_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&pfValues_device,			sizeof(pfValues_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&puKeys_device,			sizeof(puKeys_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&puiSegFlags_device,		sizeof(puiSegFlags_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&pfCoefs_device,			sizeof(pfCoefs_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&pfCoefSums_device,		sizeof(pfCoefSums_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&pfCompactedCoefs_device,	sizeof(pfCompactedCoefs_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&puCompactedKeys_device,	sizeof(puCompactedKeys_device[0]) * uMaxNrOfElementsOnTheDevice));

		CUDA_SAFE_CALL(cudaMalloc((void**)&puNrOfCompactedCoefs_device,	sizeof(puNrOfCompactedCoefs_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&puNrOfCompactedKeys_device,	sizeof(puNrOfCompactedKeys_device[0]) * uMaxNrOfElementsOnTheDevice));

		bIsInitialized = true;
	}

	void
	CCudaDWT::
	_InitEncoder
	(
		size_t				uNrOfElements,
		const uint4			pu4BinSubs[],
		const float			pfValues[],
		void* _Reserved
	)
	{
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
	}

	void
	CCudaDWT::
	_Encode
	(
		size_t				uNrOfElements,
		size_t				uNrOfDims,
		const unsigned int	puLevels[],
		const unsigned int	puWaveletLengths[],

		size_t				*puNrOfElements,
		unsigned int		puKeys_host[],
		float				pfCoefs_host[],
		void* _Reserved
	)
	{
		LIBCLOCK_INIT(bIsPrintingTiming, __FUNCTION__);

		// copy the lengths of the local coefficient array, wavelet lengths, and levels
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		CUDA_SAFE_CALL(
			cudaMemcpyToSymbol(
				"puLevels_const", 
				&puLevels[0], 
				sizeof(puLevels_const[0]) * uNrOfDims,
				0,
				cudaMemcpyHostToDevice) );

		CUDA_SAFE_CALL(
			cudaMemcpyToSymbol(
				"puWaveletLengths_const", 
				&puWaveletLengths[0], 
				sizeof(puWaveletLengths_const[0]) * uNrOfDims, 
				0,
				cudaMemcpyHostToDevice) );
		LIBCLOCK_END(bIsPrintingTiming);

		LIBCLOCK_BEGIN(bIsPrintingTiming);
		// 
		dim3 v3Blk = dim3(BLOCK_SIZE);
		dim3 v3Grid = dim3((size_t)ceilf((float)uNrOfElements / (float)v3Blk.x));

		_ProjToWavelet_kernel<<<v3Grid, v3Blk, 0>>>(
			&pu4BinSub_device[0],	// the tuples of <bin, data_subscripts> of all elements
			&pfValues_device[0],	// the counts of all elements
			(unsigned int)uNrOfDims, 
			(unsigned int)uNrOfElements,
			&puKeys_device[0],		// output: the keys of all elements. The keys are composed of bin and local_subscripts
			&pfCoefs_device[0]		// output: the wavelet projection of the current wavelet 
			);
		CUT_CHECK_ERROR("_ProjToWavelet_kernel() failed");
		LIBCLOCK_END(bIsPrintingTiming);

		// sort the wavelet projection according to the key composed by the bin and local subscripts
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		ASSERT_CUDPP(cudppSort(
			planSort,				
			&puKeys_device[0],		
			&pfCoefs_device[0],
			uNrOfElements));
		LIBCLOCK_END(bIsPrintingTiming);

		// mark the segments. the beginning of a segment is marked as 1, and all other elements are marked as 0
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		_MarkSegments_kernel<<<v3Grid, v3Blk, 0>>>(
			&puKeys_device[0],
			(unsigned int)uNrOfElements,
			&puiSegFlags_device[0]);
		CUT_CHECK_ERROR("_MarkSegments_kernel() failed");
		LIBCLOCK_END(bIsPrintingTiming);

		// compute the sum of the segments
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		ASSERT_CUDPP(cudppSegmentedScan(
			planSegScanCoefs,
			&pfCoefSums_device[0],
			&pfCoefs_device[0],
			&puiSegFlags_device[0],
			uNrOfElements));
		LIBCLOCK_END(bIsPrintingTiming);

		// compact the result
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		size_t uNrOfCompactedCoefs_host;
		ASSERT_CUDPP(cudppCompact(
			planCompactCoefs,
			&pfCompactedCoefs_device[0],
			puNrOfCompactedCoefs_device,
			&pfCoefSums_device[0],
			&puiSegFlags_device[0],
			uNrOfElements));
		LIBCLOCK_END(bIsPrintingTiming);

		LIBCLOCK_BEGIN(bIsPrintingTiming);
		CUDA_SAFE_CALL(
			cudaMemcpy(
				&uNrOfCompactedCoefs_host, 
				puNrOfCompactedCoefs_device, 
				sizeof(uNrOfCompactedCoefs_host), 
				cudaMemcpyDeviceToHost));
		LIBCLOCK_END(bIsPrintingTiming);

		LIBCLOCK_BEGIN(bIsPrintingTiming);
		CUDA_SAFE_CALL(
			cudaMemcpy(
				&pfCoefs_host[0],
				&pfCompactedCoefs_device[0], 
				uNrOfCompactedCoefs_host * sizeof(pfCoefs_host[0]),
				cudaMemcpyDeviceToHost) );
		LIBCLOCK_END(bIsPrintingTiming);

		// compact the keys
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		size_t uNrOfCompactedKeys_host;
		ASSERT_CUDPP(cudppCompact(
			planCompactKeys,
			&puCompactedKeys_device[0],
			puNrOfCompactedKeys_device,
			&puKeys_device[0],
			&puiSegFlags_device[0],
			uNrOfElements));
		LIBCLOCK_END(bIsPrintingTiming);

		LIBCLOCK_BEGIN(bIsPrintingTiming);
		CUDA_SAFE_CALL(
			cudaMemcpy(
				&uNrOfCompactedKeys_host, 
				puNrOfCompactedKeys_device, 
				sizeof(uNrOfCompactedKeys_host), 
				cudaMemcpyDeviceToHost));
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

