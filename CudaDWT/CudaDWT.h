#pragma once

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
		DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE = 0x0400000,
		#endif	// #if	defined(_DEBUG)	

		PARAMETER_END
	};

	class CCudaDWT
	{
	protected:
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

		bool bIsInitialized;
		bool bIsPrintingTiming;
	public:
		CCudaDWT():
			bIsPrintingTiming(false),
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
			unsigned int		puKeys_host[],
			float				pfCoefs_host[],
			void* _Reserved = NULL
		);

		void
		_Free
		(
			void* _Reserved = NULL
		);
	};
};

