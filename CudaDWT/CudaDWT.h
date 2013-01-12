#pragma once

// ADD-BY-LEETEN 01/11/2013-BEGIN
#define WITH_CUDA_MALLOC_HOST	0
// ADD-BY-LEETEN 01/11/2013-END

#include <cuda_runtime_api.h>
#include <cudpp.h>

#if defined(WIN32)
	#pragma comment (lib, "cudart.lib")
	#if	defined(_DEBUG)
		#if	defined(WIN64)
		#pragma comment (lib, "cudpp64d.lib")
		#endif	// #if	defined(WIN64)
		#pragma comment (lib, "CudaDWT_d.lib")
	#else	// #if	defined(_DEBUG)
		#if	defined(WIN64)
		#pragma comment (lib, "cudpp64.lib")
		#endif	// #if	defined(WIN64)
		#pragma comment (lib, "CudaDWT_r.lib")
	#endif	// #if	defined(_DEBUG)
#endif	// #if defined(WIN32)

namespace CudaDWT
{
	enum
	{
		GPU_MAX_NR_OF_DIMS = 3,

		#if	defined(_DEBUG)	
		DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE = 0x0400,
		#else	// #if	defined(_DEBUG)	
		// MOD-BY-LEETEN 01/11/2013-FROM:		DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE = 0x0400000,
		// MOD-BY-LEETEN 01/11/2013-FROM:		DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE = 0x1000000,
		DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE = 0x8000000,
		// MOD-BY-LEETEN 01/11/2013-END
		// MOD-BY-LEETEN 01/11/2013-END
		#endif	// #if	defined(_DEBUG)	

		PARAMETER_END
	};

	class CCudaDWT
	{
	protected:
		// ADD-BY-LEETEN 01/11/2013-BEGIN
		dim3 v3Blk;

		dim3 v3Grid;
		// ADD-BY-LEETEN 01/11/2013-END

		// CUDPP handles
		CUDPPHandle theCudpp;

		CUDPPConfiguration configSort;
		CUDPPHandle planSort;

		CUDPPConfiguration configSegScanCoefs;
		CUDPPHandle planSegScanCoefs;

		size_t *puNrOfCompactedKeys_device;
		CUDPPConfiguration configCompactCoefs;
		CUDPPHandle planCompactCoefs;

		size_t *puNrOfCompactedCoefs_device;
		CUDPPConfiguration configCompactKeys;
		CUDPPHandle planCompactKeys;

		//
		uint4 *pu4BinSub_device;
		float* pfValues_device;
		unsigned int* puKeys_device;
		unsigned int* puiSegFlags_device;
		float* pfCoefs_device;
		float* pfCoefSums_device;
		float* pfCompactedCoefs_device;
		unsigned int* puCompactedKeys_device;

		// ADD-BY-LEETEN 01/11/2013-BEGIN
		#if	WITH_CUDA_MALLOC_HOST 
		size_t *puNrOfCompactedKeys_host;
		size_t *puNrOfCompactedCoefs_host;
		unsigned int* puKeys_host;
		float*	pfCoefs_host;
		#endif	// #if	WITH_CUDA_MALLOC_HOST 
		// ADD-BY-LEETEN 01/11/2013-END

		bool bIsInitialized;
		// DEL-BY-LEETEN 01/11/2013:	bool bIsPrintingTiming;
	public:
		CCudaDWT():
			// DEL-BY-LEETEN 01/11/2013:	bIsPrintingTiming(false),
			bIsInitialized(false)
		{
		}

		virtual 
			~CCudaDWT();

		void
		_Init
		(
			size_t* puMaxNrOfElementsOnTheDevice,
			void* _Reserved = NULL
		);

		void
		_InitEncoder
		(
			size_t				uNrOfElements,
			const uint4			pu4BinSubs[],
			const float			pfValues[],
			void* _Reserved = NULL
		);

		void
		_Encode
		(
			size_t				uNrOfElements,
			size_t				uNrOfDims,
			const unsigned int	puLevels[],
			const unsigned int	puWaveletLengths[],

			size_t				*puNrOfElements,
			#if	!WITH_CUDA_MALLOC_HOST	// ADD-BY-LEETEN 01/11/2013
			unsigned int		puKeys_host[],
			float				pfCoefs_host[],
			// ADD-BY-LEETEN 01/11/2013-BEGIN
			#else	// #if	!WITH_CUDA_MALLOC_HOST
			unsigned int		puKeys[],
			float				pfCoefs[],
			#endif	// #if	!WITH_CUDA_MALLOC_HOST
			// ADD-BY-LEETEN 01/11/2013-END
			int iTimingPrintingLevel = 0,	// ADD-BY-LEETEN 01/11/2013
			void* _Reserved = NULL
		);

		void
		_Free
		(
			void* _Reserved = NULL
		);
	};
};

