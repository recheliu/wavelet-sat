#include <vector>
using namespace std;
#include <stdio.h>
#include <assert.h>
#include "libclock.h"
#include "liblog.h"
#include "cuda_macro.h"
#include "CudaDWT.h"

// ADD-BY-LEETEN 03/29/2013-BEGIN
// DEL-BY-LEETEN 04/08/2013:	#if	WITH_DOUBLE_COEF
#if	0	// DEL-BY-LEETEN 04/08/2013-BEGIN
#if	WITH_DOUBLE_COEF
typedef double typeCoef;
#else
typedef float typeCoef;
#endif	// #if	WITH_DOUBLE_COEF
#endif	// DEL-BY-LEETEN 04/08/2013-END
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
// DEL-BY-LEETEN 04/08/2013:	#endif	// #if	WITH_DOUBLE_COEF	
// ADD-BY-LEETEN 03/29/2013-END

#if	0	// DEL-BY-LEETEN 04/08/2013-BEGIN
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

#endif	// DEL-BY-LEETEN 04/08/2013-END

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
			// DEL-BY-LEETEN 04/08/2013:	FREE_MEMORY(puiSegFlags_device);
			FREE_MEMORY(pfCoefs_device);
			// DEL-BY-LEETEN 04/08/2013:	FREE_MEMORY(pfCoefSums_device);
			FREE_MEMORY(pfCompactedCoefs_device);
			FREE_MEMORY(puCompactedKeys_device);

			// ADD-BY-LEETEN 01/13/2013-BEGIN
			FREE_MEMORY(puOnes_device);
			// DEL-BY-LEETEN 04/08/2013:	FREE_MEMORY(puSegCounts_device);
			FREE_MEMORY(puCompactedSegCounts_device);
			// DEL-BY-LEETEN 04/08/2013:	FREE_MEMORY(puNrOfCompactedSegCounts_device);
			// ADD-BY-LEETEN 01/13/2013-END

			// ADD-BY-LEETEN 01/11/2013-BEGIN
			#if		WITH_CUDA_MALLOC_HOST	
			FREE_MEMORY_ON_HOST(puNrOfCompactedKeys_host);
			FREE_MEMORY_ON_HOST(puNrOfCompactedCoefs_host);
			FREE_MEMORY_ON_HOST(pfCoefs_host);
			FREE_MEMORY_ON_HOST(puKeys_host);
			#endif	// #if		WITH_CUDA_MALLOC_HOST	
			// ADD-BY-LEETEN 01/11/2013-END

			#if	0	// DEL-BY-LEETEN 04/08/2013-BEGIN
			FREE_MEMORY(puNrOfCompactedCoefs_device);
			FREE_MEMORY(puNrOfCompactedKeys_device);

			// free cudpp resources
			#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
			if( planSort )
				ASSERT_CUDPP(cudppDestroyPlan(planSort));  
			#endif	// #if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
			if( planSegScanCoefs )
				ASSERT_CUDPP(cudppDestroyPlan(planSegScanCoefs));
			if( planCompactCoefs )
				ASSERT_CUDPP(cudppDestroyPlan(planCompactCoefs));
			if( planCompactKeys )
				ASSERT_CUDPP(cudppDestroyPlan(planCompactKeys));
			// ADD-BY-LEETEN 01/13/2013-BEGIN
			if( planSegScanCounts )
				ASSERT_CUDPP(cudppDestroyPlan(planSegScanCounts));

			if( planCompactSegCounts )
				ASSERT_CUDPP(cudppDestroyPlan(planCompactSegCounts));
			// ADD-BY-LEETEN 01/13/2013-END
			if( theCudpp )
				ASSERT_CUDPP(cudppDestroy(theCudpp));
			#endif	// DEL-BY-LEETEN 04/08/2013-END
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
		// DEL-BY-LEETEN 04/08/2013:	ASSERT_CUDPP(cudppCreate(&theCudpp));
		if( *puMaxNrOfElementsOnTheDevice > CudaDWT::DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE )
		{
			LOG_ERROR(cerr<<"uMaxNrOfElementsOnTheDevice is clampped to CudaDWT::DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE");
			*puMaxNrOfElementsOnTheDevice = CudaDWT::DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE;
		}

		size_t uMaxNrOfElementsOnTheDevice = *puMaxNrOfElementsOnTheDevice;
		#if	0	// DEL-BY-LEETEN 04/08/2013-BEGIN
		#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
		configSort.op = CUDPP_ADD;
		configSort.datatype = CUDPP_UINT;
		configSort.algorithm = CUDPP_SORT_RADIX;
		configSort.options = CUDPP_OPTION_KEY_VALUE_PAIRS;
		ASSERT_CUDPP(cudppPlan(theCudpp, &planSort, configSort, uMaxNrOfElementsOnTheDevice, 1, 0));  
		#endif	// #if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
		configSegScanCoefs.op = CUDPP_ADD;
		#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
		configSegScanCoefs.datatype = CUDPP_FLOAT;
		// ADD-BY-LEETEN 03/29/2013-BEGIN
		#else	// #if	!WITH_DOUBLE_COEF
		configSegScanCoefs.datatype = CUDPP_DOUBLE;
		#endif	// #if	!WITH_DOUBLE_COEF
		// ADD-BY-LEETEN 03/29/2013-END
		configSegScanCoefs.algorithm = CUDPP_SEGMENTED_SCAN;
		configSegScanCoefs.options = CUDPP_OPTION_BACKWARD | CUDPP_OPTION_INCLUSIVE;
		ASSERT_CUDPP(cudppPlan(theCudpp, &planSegScanCoefs, configSegScanCoefs, uMaxNrOfElementsOnTheDevice, 1, 0));

		#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
		configCompactCoefs.datatype = CUDPP_FLOAT;
		// ADD-BY-LEETEN 03/29/2013-BEGIN
		#else	// #if	!WITH_DOUBLE_COEF	
		configCompactCoefs.datatype = CUDPP_DOUBLE;
		#endif	// #if	!WITH_DOUBLE_COEF
		// ADD-BY-LEETEN 03/29/2013-END
		configCompactCoefs.algorithm = CUDPP_COMPACT;
		configCompactCoefs.options = CUDPP_OPTION_FORWARD;
		ASSERT_CUDPP(cudppPlan(theCudpp, &planCompactCoefs, configCompactCoefs, uMaxNrOfElementsOnTheDevice, 1, 0));

		configCompactKeys.datatype = CUDPP_UINT;
		configCompactKeys.algorithm = CUDPP_COMPACT;
		configCompactKeys.options = CUDPP_OPTION_FORWARD;
		ASSERT_CUDPP(cudppPlan(theCudpp, &planCompactKeys, configCompactKeys, uMaxNrOfElementsOnTheDevice, 1, 0));
		#endif	// DEL-BY-LEETEN 04/08/2013-END
		// allocate the memory space
		CUDA_SAFE_CALL(cudaMalloc((void**)&pu4BinSub_device,		sizeof(pu4BinSub_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&pfValues_device,			sizeof(pfValues_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&puKeys_device,			sizeof(puKeys_device[0]) * uMaxNrOfElementsOnTheDevice));
		// DEL-BY-LEETEN 04/08/2013:	CUDA_SAFE_CALL(cudaMalloc((void**)&puiSegFlags_device,		sizeof(puiSegFlags_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&pfCoefs_device,			sizeof(pfCoefs_device[0]) * uMaxNrOfElementsOnTheDevice));
		// DEL-BY-LEETEN 04/08/2013:	CUDA_SAFE_CALL(cudaMalloc((void**)&pfCoefSums_device,		sizeof(pfCoefSums_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&pfCompactedCoefs_device,	sizeof(pfCompactedCoefs_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&puCompactedKeys_device,	sizeof(puCompactedKeys_device[0]) * uMaxNrOfElementsOnTheDevice));

		#if	0	// DEL-BY-LEETEN 04/08/2013-BEGIN
		// ADD-BY-LEETEN 01/13/2013-BEGIN
		configSegScanCounts.datatype = CUDPP_UINT;
		configSegScanCounts.algorithm = CUDPP_SEGMENTED_SCAN;
		configSegScanCounts.options = CUDPP_OPTION_BACKWARD | CUDPP_OPTION_INCLUSIVE;
		ASSERT_CUDPP(cudppPlan(theCudpp, &planSegScanCounts, configSegScanCounts, uMaxNrOfElementsOnTheDevice, 1, 0));

		configCompactSegCounts.datatype = CUDPP_UINT;
		configCompactSegCounts.algorithm = CUDPP_COMPACT;
		configCompactSegCounts.options = CUDPP_OPTION_FORWARD;
		ASSERT_CUDPP(cudppPlan(theCudpp, &planCompactSegCounts, configCompactSegCounts, uMaxNrOfElementsOnTheDevice, 1, 0));
		#endif	// DEL-BY-LEETEN 04/08/2013-END
		CUDA_SAFE_CALL(cudaMalloc((void**)&puOnes_device,			sizeof(puOnes_device[0]) * uMaxNrOfElementsOnTheDevice));
		// DEL-BY-LEETEN 04/08/2013:	CUDA_SAFE_CALL(cudaMalloc((void**)&puSegCounts_device,		sizeof(puSegCounts_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&puCompactedSegCounts_device,			sizeof(puCompactedSegCounts_device[0]) * uMaxNrOfElementsOnTheDevice));
		// DEL-BY-LEETEN 04/08/2013:	CUDA_SAFE_CALL(cudaMalloc((void**)&puNrOfCompactedSegCounts_device,		sizeof(puNrOfCompactedSegCounts_device[0])));
		vector<unsigned int> vuOnes;
		vuOnes.assign(uMaxNrOfElementsOnTheDevice, 1);
		CUDA_SAFE_CALL(cudaMemcpy(puOnes_device, vuOnes.data(), sizeof(puOnes_device[0]) * vuOnes.size(), cudaMemcpyHostToDevice));
		// ADD-BY-LEETEN 01/13/2013-END

		// ADD-BY-LEETEN 01/11/2013-BEGIN
		#if		WITH_CUDA_MALLOC_HOST	
		CUDA_SAFE_CALL(cudaMallocHost((void**)&puNrOfCompactedKeys_host,			sizeof(puNrOfCompactedKeys_host[0])));
		CUDA_SAFE_CALL(cudaMallocHost((void**)&puNrOfCompactedCoefs_host,			sizeof(puNrOfCompactedCoefs_host[0])));
		CUDA_SAFE_CALL(cudaMallocHost((void**)&puKeys_host,			sizeof(puKeys_host[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMallocHost((void**)&pfCoefs_host,		sizeof(pfCoefs_host[0]) * uMaxNrOfElementsOnTheDevice));
		#endif	// #if	WITH_CUDA_MALLOC_HOST	
		// ADD-BY-LEETEN 01/11/2013-END
		#if	0	// DEL-BY-LEETEN 04/08/2013-BEGIN
		CUDA_SAFE_CALL(cudaMalloc((void**)&puNrOfCompactedCoefs_device,	sizeof(puNrOfCompactedCoefs_device[0]) * uMaxNrOfElementsOnTheDevice));
		CUDA_SAFE_CALL(cudaMalloc((void**)&puNrOfCompactedKeys_device,	sizeof(puNrOfCompactedKeys_device[0]) * uMaxNrOfElementsOnTheDevice));
		#endif	// DEL-BY-LEETEN 04/08/2013-END
		bIsInitialized = true;
	}

	void
	CCudaDWT::
	_InitEncoder
	(
		// ADD-BY-LEETEN 01/18/2012-BEGIN
		size_t uNrOfDims,
		unsigned int puCoefLengths[],
		// ADD-BY-LEETEN 01/18/2012-END
		size_t				uNrOfElements,
		const uint4			pu4BinSubs[],
		#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
		const float			pfValues[],
		// ADD-BY-LEETEN 03/29/2013-BEGIN
		#else	// #if	!WITH_DOUBLE_COEF
		// MOD-BY-LEETEN 04/08/2013-FROM:		const double		pfValues[],
		const typeCoef	pfValues[],
		// MOD-BY-LEETEN 04/08/2013-END
		#endif	// #if	!WITH_DOUBLE_COEF
		// ADD-BY-LEETEN 03/29/2013-END
		bool bWithCpuBucketSort,	// ADD-BY-LEETEN 01/13/2013
		void* _Reserved
	)
	{
		this->bWithCpuBucketSort = bWithCpuBucketSort;	// ADD-BY-LEETEN 01/13/2012

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

		// ADD-BY-LEETEN 01/18/2012-BEGIN
		CUDA_SAFE_CALL(
			cudaMemcpyToSymbol(
				puCoefLengths_const, 
				puCoefLengths,
				uNrOfDims * sizeof(puCoefLengths_const[0]), 
				0,
				cudaMemcpyHostToDevice));
		// ADD-BY-LEETEN 01/18/2012-END

		// ADD-BY-LEETEN 01/11/2013-BEGIN
		v3Blk = dim3(BLOCK_SIZE);
		size_t uNrOfBlocks = (size_t)ceilf((float)uNrOfElements / (float)v3Blk.x);
		size_t uGridSizeX = (size_t)ceil(sqrtf((float)uNrOfBlocks));
		size_t uGridSizeY = (size_t)ceil((float)uNrOfBlocks/(float)uGridSizeX);
		v3Grid = dim3( (unsigned int)uGridSizeX, (unsigned int)uGridSizeY );
		// ADD-BY-LEETEN 01/11/2013-END
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
		#if		!WITH_CUDA_MALLOC_HOST	// ADD-BY-LEETEN 01/11/2013
		unsigned int		puKeys_host[],
		#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
		float				pfCoefs_host[],
		// ADD-BY-LEETEN 03/29/2013-BEGIN
		#else	// #if	!WITH_DOUBLE_COEF
		// MOD-BY-LEETEN 04/08/2013-FROM:		double				pfCoefs_host[],
		typeCoef			pfCoefs_host[],
		// MOD-BY-LEETEN 04/08/2013-END
		#endif	// #if	!WITH_DOUBLE_COEF
		// ADD-BY-LEETEN 03/29/2013-END
		// ADD-BY-LEETEN 01/11/2013-BEGIN
		#else	// #if		!WITH_CUDA_MALLOC_HOST
		unsigned int		puKeys[],
		float				pfCoefs[],
		#endif	// #if		!WITH_CUDA_MALLOC_HOST
		// ADD-BY-LEETEN 01/11/2013-END

		unsigned int		puSegCounts_host[],	// ADD-BY-LEETEN 01/13/2013

		int iTimingPrintingLevel,	// ADD-BY-LEETEN 01/11/2013
		void* _Reserved
	)
	{
		bool bIsPrintingTiming = (iTimingPrintingLevel > 0)?true:false;	// ADD-BY-LEETEN 01/11/2013
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
			(unsigned int)uNrOfDims, 
			(unsigned int)uNrOfElements,
			&puKeys_device[0],		// output: the keys of all elements. The keys are composed of bin and local_subscripts
			&pfCoefs_device[0]		// output: the wavelet projection of the current wavelet 
			);
		CUT_CHECK_ERROR("_ProjToWavelet_kernel() failed");
		LIBCLOCK_END(bIsPrintingTiming);

		// ADD-BY-LEETEN 01/13/2013-BEGIN
		if( !bWithCpuBucketSort )	
		{
		// ADD-BY-LEETEN 01/13/2013-END
		// sort the wavelet projection according to the key composed by the bin and local subscripts
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		#if	0	// MOD-BY-LEETEN 04/08/2013-FROM:
		#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
		ASSERT_CUDPP(cudppSort(
			planSort,				
			&puKeys_device[0],		
			&pfCoefs_device[0],
			uNrOfElements));
		// ADD-BY-LEETEN 03/29/2013-BEGIN
		#else	// #if	!WITH_DOUBLE_COEF
		thrust::device_ptr<unsigned int> vKeys_device(puKeys_device);
		thrust::device_ptr<double> vCoefs_device(pfCoefs_device);
		thrust::sort_by_key(vKeys_device, vKeys_device + uNrOfElements, vCoefs_device);
		#endif	// #if	!WITH_DOUBLE_COEF
		#else	// MOD-BY-LEETEN 04/08/2013-TO:
		thrust::device_ptr<unsigned int> vKeys_device(puKeys_device);
		thrust::device_ptr<typeCoef> vCoefs_device(pfCoefs_device);
		thrust::sort_by_key(vKeys_device, vKeys_device + uNrOfElements, vCoefs_device);
		#endif	// MOD-BY-LEETEN 04/08/2013-END
		// ADD-BY-LEETEN 03/29/2013-END
		LIBCLOCK_END(bIsPrintingTiming);
		}	// ADD-BY-LEETEN 01/13/2013

		#if	0	// MOD-BY-LEETEN 04/08/2013-FROM:
		// mark the segments. the beginning of a segment is marked as 1, and all other elements are marked as 0
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		_MarkSegments_kernel<<<v3Grid, v3Blk, 0>>>(
			&puKeys_device[0],
			(unsigned int)uNrOfElements,
			&puiSegFlags_device[0]);
		CUT_CHECK_ERROR("_MarkSegments_kernel() failed");
		LIBCLOCK_END(bIsPrintingTiming);

		// ADD-BY-LEETEN 01/13/2013-BEGIN
		// compute the count per segment
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		ASSERT_CUDPP(cudppSegmentedScan(
			planSegScanCounts,
			&puSegCounts_device[0],
			&puOnes_device[0],
			&puiSegFlags_device[0],
			uNrOfElements));
		LIBCLOCK_END(bIsPrintingTiming);

		// compact the result
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		ASSERT_CUDPP(cudppCompact(
			planCompactSegCounts,
			&puCompactedSegCounts_device[0],
			puNrOfCompactedSegCounts_device,
			&puSegCounts_device[0],
			&puiSegFlags_device[0],
			uNrOfElements));
		LIBCLOCK_END(bIsPrintingTiming);

		LIBCLOCK_BEGIN(bIsPrintingTiming);
		size_t uNrOfCompactedSegCounts_host;
		CUDA_SAFE_CALL(
			cudaMemcpy(
				&uNrOfCompactedSegCounts_host, 
				puNrOfCompactedSegCounts_device, 
				sizeof(uNrOfCompactedSegCounts_host), 
				cudaMemcpyDeviceToHost));
		LIBCLOCK_END(bIsPrintingTiming);
		#else	// MOD-BY-LEETEN 04/08/2013-TO:
		thrust::device_ptr<unsigned int> vuKeys_device(puKeys_device);
		thrust::device_ptr<unsigned int> vuOnes_device(puOnes_device);
		thrust::device_ptr<unsigned int> vuCompactedKeys_device(puCompactedKeys_device);
		thrust::device_ptr<unsigned int> vuCompactedSegCounts_device(puCompactedSegCounts_device);
		thrust::pair< thrust::device_ptr<unsigned int>, thrust::device_ptr<unsigned int> > pairEnds = 
			thrust::reduce_by_key<
				thrust::device_ptr<unsigned int>, 
				thrust::device_ptr<unsigned int>, 
				thrust::device_ptr<unsigned int>, 
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
		#endif	// MOD-BY-LEETEN 04/08/2013-END
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		CUDA_SAFE_CALL(
			cudaMemcpy(
				&puSegCounts_host[0],
				&puCompactedSegCounts_device[0], 
				uNrOfCompactedSegCounts_host * sizeof(puSegCounts_host[0]),
				cudaMemcpyDeviceToHost) );
		LIBCLOCK_END(bIsPrintingTiming);
		// ADD-BY-LEETEN 01/13/2013-END
		#if	0	// MOD-BY-LEETEN 04/08/2013-FROM:
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
		#if		!WITH_CUDA_MALLOC_HOST	// ADD-BY-LEETEN 01/11/2013
		CUDA_SAFE_CALL(
			cudaMemcpy(
				&uNrOfCompactedCoefs_host, 
				puNrOfCompactedCoefs_device, 
				sizeof(uNrOfCompactedCoefs_host), 
				cudaMemcpyDeviceToHost));
		// ADD-BY-LEETEN 01/11/2013-BEGIN
		#else	// #if		!WITH_CUDA_MALLOC_HOST	
		CUDA_SAFE_CALL(
			cudaMemcpy(
				puNrOfCompactedCoefs_host, 
				puNrOfCompactedCoefs_device, 
				sizeof(uNrOfCompactedCoefs_host), 
				cudaMemcpyDeviceToHost));
		uNrOfCompactedCoefs_host = *puNrOfCompactedCoefs_host;
		#endif	// #if		!WITH_CUDA_MALLOC_HOST	
		// ADD-BY-LEETEN 01/11/2013-END
		LIBCLOCK_END(bIsPrintingTiming);
		#else	// MOD-BY-LEETEN 04/08/2013-TO::
		thrust::device_ptr<typeCoef> vfCoefs_device(pfCoefs_device);
		thrust::device_ptr<typeCoef> vfCompactedCoefs_device(pfCompactedCoefs_device);
		size_t uNrOfCompactedCoefs_host;
		{
			thrust::pair<thrust::device_ptr<unsigned int>, thrust::device_ptr<typeCoef> > pairEnds = 
				thrust::reduce_by_key<
					thrust::device_ptr<unsigned int>, 
					thrust::device_ptr<typeCoef>, 
					thrust::device_ptr<unsigned int>, 
					thrust::device_ptr<typeCoef> >
				(
					vuKeys_device, 
					vuKeys_device + uNrOfElements, 
					vfCoefs_device,
					vuCompactedKeys_device,
					vfCompactedCoefs_device);
			uNrOfCompactedCoefs_host = size_t(pairEnds.second - vfCompactedCoefs_device);
		}
		#endif	// MOD-BY-LEETEN 04/08/2013-END
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		CUDA_SAFE_CALL(
			cudaMemcpy(
				&pfCoefs_host[0],
				&pfCompactedCoefs_device[0], 
				uNrOfCompactedCoefs_host * sizeof(pfCoefs_host[0]),
				cudaMemcpyDeviceToHost) );
		// ADD-BY-LEETEN 01/11/2013-BEGIN
		#if		WITH_CUDA_MALLOC_HOST	
		memcpy(&pfCoefs[0], &pfCoefs_host[0], uNrOfCompactedCoefs_host * sizeof(pfCoefs[0]));
		#endif	// #if		WITH_CUDA_MALLOC_HOST	
		// ADD-BY-LEETEN 01/11/2013-END
		LIBCLOCK_END(bIsPrintingTiming);

		#if	0	// DEL-BY-LEETEN 04/08/2013-BEGIN
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
		#if		!WITH_CUDA_MALLOC_HOST	// ADD-BY-LEETEN 01/11/2013
		CUDA_SAFE_CALL(
			cudaMemcpy(
				&uNrOfCompactedKeys_host, 
				puNrOfCompactedKeys_device, 
				sizeof(uNrOfCompactedKeys_host), 
				cudaMemcpyDeviceToHost));
		// ADD-BY-LEETEN 01/11/2013-BEGIN
		#else	// #if		!WITH_CUDA_MALLOC_HOST	
		CUDA_SAFE_CALL(
			cudaMemcpy(
				puNrOfCompactedKeys_host, 
				puNrOfCompactedKeys_device, 
				sizeof(uNrOfCompactedKeys_host), 
				cudaMemcpyDeviceToHost));
		uNrOfCompactedKeys_host = *puNrOfCompactedKeys_host;
		#endif	// #if		!WITH_CUDA_MALLOC_HOST	
		// ADD-BY-LEETEN 01/11/2013-END
		LIBCLOCK_END(bIsPrintingTiming);
		#endif	// DEL-BY-LEETEN 04/08/2013-END
		// download the keys and the coefficinets back 
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		CUDA_SAFE_CALL(
			cudaMemcpy(
				&puKeys_host[0],
				&puCompactedKeys_device[0], 
				uNrOfCompactedKeys_host * sizeof(puKeys_host[0]),
				cudaMemcpyDeviceToHost) );
		// ADD-BY-LEETEN 01/11/2013-BEGIN
		#if		WITH_CUDA_MALLOC_HOST	
		memcpy(&puKeys[0], &puKeys_host[0], uNrOfCompactedKeys_host * sizeof(puKeys[0]));
		#endif	// #if		WITH_CUDA_MALLOC_HOST	
		// ADD-BY-LEETEN 01/11/2013-END
		LIBCLOCK_END(bIsPrintingTiming);

		ASSERT_OR_LOG(uNrOfCompactedKeys_host == uNrOfCompactedCoefs_host, cerr<<"Unmatched #keys and #coefs.");
		*puNrOfElements = uNrOfCompactedKeys_host;
		LIBCLOCK_PRINT(bIsPrintingTiming);
	}
};

