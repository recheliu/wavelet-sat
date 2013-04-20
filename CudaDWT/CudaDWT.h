#pragma once

#include <cuda_runtime_api.h>

#if defined(WIN32)
	#pragma comment (lib, "cudart.lib")
	#if	defined(_DEBUG)
		#pragma comment (lib, "CudaDWT_d.lib")
	#else	// #if	defined(_DEBUG)
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

		//
		uint4 *pu4BinSub_device;
		#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
		float* pfValues_device;
		// ADD-BY-LEETEN 03/29/2013-BEGIN
		#else	// #if	!WITH_DOUBLE_COEF	
		typeCoef *pfValues_device;
		#endif	// #if	!WITH_DOUBLE_COEF	
		// ADD-BY-LEETEN 03/29/2013-END
		unsigned int* puKeys_device;

		typeCoef* pfCoefs_device;
		typeCoef* pfCompactedCoefs_device;
		unsigned int* puCompactedKeys_device;

		unsigned int *puOnes_device;
		unsigned int *puCompactedSegCounts_device;
		// ADD-BY-LEETEN 01/13/2013-END

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
			unsigned int		puKeys_host[],
			#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
			float				pfCoefs_host[],
			// ADD-BY-LEETEN 03/29/2013-BEGIN
			#else	// #if	!WITH_DOUBLE_COEF
			double				pfCoefs_host[],
			#endif	// #if	!WITH_DOUBLE_COEF
			// ADD-BY-LEETEN 03/29/2013-END

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

