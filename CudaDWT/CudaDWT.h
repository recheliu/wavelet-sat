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


#include "Base.h"

namespace CudaDWT
{
	typedef WaveletSAT::typeWavelet	typeCoef;
	typedef	double	typeValue;		
	typedef unsigned long long typeKey;	
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
		dim3 v3Blk;

		dim3 v3Grid;

		//
		uint4 *pu4BinSub_device;
		typeValue	*pfValues_device;
		typeKey *puKeys_device;

		typeCoef* pfCoefs_device;
		typeCoef* pfCompactedCoefs_device;
		typeKey* puCompactedKeys_device;

		unsigned int *puOnes_device;
		unsigned int *puCompactedSegCounts_device;

		bool bIsInitialized;

		bool bWithCpuBucketSort;
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
			size_t uNrOfDims,
			unsigned int puCoefLengths[],
			size_t				uNrOfElements,
			const uint4			pu4BinSubs[],
			const typeValue pfValues[],
			bool bWithCpuBucketSort = false,	
			void* _Reserved = NULL
		);

		void
		_Encode
		(
			size_t				uNrOfBins,	
			size_t				uNrOfElements,
			size_t				uNrOfDims,
			const unsigned int	puLevels[],
			const unsigned int	puWaveletLengths[],

			size_t				*puNrOfElements,
			typeKey		puKeys_host[],
			typeCoef	pfCoefs_host[],

			unsigned int		puSegCounts_host[],	

			int iTimingPrintingLevel = 0,	
			void* _Reserved = NULL
		);

		void
		_Free
		(
			void* _Reserved = NULL
		);
	};
};


/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
