#pragma once

// ADD-BY-LEETEN 01/11/2013-BEGIN
#define WITH_CUDA_MALLOC_HOST	0
// ADD-BY-LEETEN 01/11/2013-END

#include <cuda_runtime_api.h>
// DEL-BY-LEETEN 04/08/2013:	#include <cudpp.h>

#if defined(WIN32)
	#pragma comment (lib, "cudart.lib")
	#if	defined(_DEBUG)
		#if	0	// DEL-BY-LEETEN 04/08/2013-BEGIN
		#if	defined(WIN64)
		#pragma comment (lib, "cudpp64d.lib")
		#endif	// #if	defined(WIN64)
		#endif	// DEL-BY-LEETEN 04/08/2013-END
		#pragma comment (lib, "CudaDWT_d.lib")
	#else	// #if	defined(_DEBUG)
		#if	0	// DEL-BY-LEETEN 04/08/2013-BEGIN
		#if	defined(WIN64)
		#pragma comment (lib, "cudpp64.lib")
		#endif	// #if	defined(WIN64)
		#endif	// DEL-BY-LEETEN 04/08/2013-END
		#pragma comment (lib, "CudaDWT_r.lib")
	#endif	// #if	defined(_DEBUG)
#endif	// #if defined(WIN32)

namespace CudaDWT
{
	// ADD-BY-LEETEN 04/08/2013-BEGIN
	#if	WITH_DOUBLE_COEF
	typedef double typeCoef;
	#else
	typedef float typeCoef;
	#endif	// #if	WITH_DOUBLE_COEF
	// ADD-BY-LEETEN 04/08/2013-END
	enum
	{
		GPU_MAX_NR_OF_DIMS = 3,

		#if	defined(_DEBUG)	
		DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE = 0x0400,
		#else	// #if	defined(_DEBUG)	
		DEFAULT_MAX_NR_OF_ELEMENTS_ON_THE_DEVICE = 0x8000000,
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

		#if	0	// DEL-BY-LEETEN 04/08/2013-BEGIN
		// CUDPP handles
		CUDPPHandle theCudpp;

		#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
		CUDPPConfiguration configSort;
		CUDPPHandle planSort;
		#endif	// #if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013

		CUDPPConfiguration configSegScanCoefs;
		CUDPPHandle planSegScanCoefs;

		size_t *puNrOfCompactedKeys_device;
		CUDPPConfiguration configCompactCoefs;
		CUDPPHandle planCompactCoefs;

		size_t *puNrOfCompactedCoefs_device;
		CUDPPConfiguration configCompactKeys;
		CUDPPHandle planCompactKeys;
		#endif	// DEL-BY-LEETEN 04/08/2013-END

		//
		uint4 *pu4BinSub_device;
		#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
		float* pfValues_device;
		// ADD-BY-LEETEN 03/29/2013-BEGIN
		#else	// #if	!WITH_DOUBLE_COEF	
		// MOD-BY-LEETEN 04/09/2013-FROM:		double* pfValues_device;
		typeCoef *pfValues_device;
		// MOD-BY-LEETEN 04/09/2013-END
		#endif	// #if	!WITH_DOUBLE_COEF	
		// ADD-BY-LEETEN 03/29/2013-END
		unsigned int* puKeys_device;
		// DEL-BY-LEETEN 04/08/2013:		unsigned int* puiSegFlags_device;

		#if	0	// MOD-BY-LEETEN 04/08/2013-FROM:
		#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
		float* pfCoefs_device;
		float* pfCoefSums_device;
		float* pfCompactedCoefs_device;
		// ADD-BY-LEETEN 03/29/2013-BEGIN
		#else	// #if	!WITH_DOUBLE_COEF
		double* pfCoefs_device;
		double* pfCoefSums_device;
		double* pfCompactedCoefs_device;
		#endif	// #if	!WITH_DOUBLE_COEF
		// ADD-BY-LEETEN 03/29/2013-END
		#else	// MOD-BY-LEETEN 04/08/2013-TO:
		typeCoef* pfCoefs_device;
		typeCoef* pfCompactedCoefs_device;
		#endif	// MOD-BY-LEETEN 04/08/2013-END
		unsigned int* puCompactedKeys_device;

		#if	0	// DEL-BY-LEETEN 04/08/2013-BEGIN
		// ADD-BY-LEETEN 01/13/2013-BEGIN
		CUDPPConfiguration configSegScanCounts;
		CUDPPHandle planSegScanCounts;

		size_t *puNrOfCompactedSegCounts_device;
		CUDPPConfiguration configCompactSegCounts;
		CUDPPHandle planCompactSegCounts;
		#endif	// DEL-BY-LEETEN 04/08/2013-END

		unsigned int *puOnes_device;
		// DEL-BY-LEETEN 04/08/2013:	unsigned int *puSegCounts_device;
		unsigned int *puCompactedSegCounts_device;
		// ADD-BY-LEETEN 01/13/2013-END

		// ADD-BY-LEETEN 01/11/2013-BEGIN
		#if	WITH_CUDA_MALLOC_HOST 
		size_t *puNrOfCompactedKeys_host;
		size_t *puNrOfCompactedCoefs_host;
		unsigned int* puKeys_host;
		float*	pfCoefs_host;
		#endif	// #if	WITH_CUDA_MALLOC_HOST 
		// ADD-BY-LEETEN 01/11/2013-END

		bool bIsInitialized;

		bool bWithCpuBucketSort;	// ADD-BY-LEETEN 01/13/2012
	public:
		CCudaDWT():
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
			// ADD-BY-LEETEN 01/18/2012-BEGIN
			size_t uNrOfDims,
			unsigned int puCoefLengths[],
			// ADD-BY-LEETEN 01/18/2012-END
			size_t				uNrOfElements,
			const uint4			pu4BinSubs[],
			#if	!WITH_DOUBLE_COEF				// ADD-BY-LEETEN 03/29/2013
			const float			pfValues[],
			// ADD-BY-LEETEN 03/29/2013-BEGIN
			#else	// #if	!WITH_DOUBLE_COEF
			const double			pfValues[],
			#endif	// #if	!WITH_DOUBLE_COEF
			// ADD-BY-LEETEN 03/29/2013-END
			bool bWithCpuBucketSort = false,	// ADD-BY-LEETEN 01/13/2013
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
			#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
			float				pfCoefs_host[],
			// ADD-BY-LEETEN 03/29/2013-BEGIN
			#else	// #if	!WITH_DOUBLE_COEF
			double				pfCoefs_host[],
			#endif	// #if	!WITH_DOUBLE_COEF
			// ADD-BY-LEETEN 03/29/2013-END
			// ADD-BY-LEETEN 01/11/2013-BEGIN
			#else	// #if	!WITH_CUDA_MALLOC_HOST
			unsigned int		puKeys[],
			float				pfCoefs[],
			#endif	// #if	!WITH_CUDA_MALLOC_HOST
			// ADD-BY-LEETEN 01/11/2013-END

			unsigned int		puSegCounts_host[],	// ADD-BY-LEETEN 01/13/2013

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

