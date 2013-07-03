__constant__ unsigned int puWaveletLengths_const[CudaDWT::GPU_MAX_NR_OF_DIMS];
__constant__ unsigned int puLevels_const[CudaDWT::GPU_MAX_NR_OF_DIMS];
__constant__ unsigned int puCoefLengths_const[CudaDWT::GPU_MAX_NR_OF_DIMS];	// ADD-BY-LEETEN 01/18/2012

__global__ 
void 
_ProjToWavelet_kernel
(
	const uint4 pu4BinSub_device[],	// the tuples of <bin, data_subscripts> of all elements
	#if	0	// MOD-BY-LEETEN 2013/07/03-FROM:
	#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
	const float pfCounts_device[],	// the counts of all elements
	// ADD-BY-LEETEN 03/29/2013-BEGIN
	#else	// #if	!WITH_DOUBLE_COEF
	const double pfCounts_device[],	// the counts of all elements
	#endif	// #if	!WITH_DOUBLE_COEF
	// ADD-BY-LEETEN 03/29/2013-END
	#else	// MOD-BY-LEETEN 2013/07/03-TO:
	const CudaDWT::typeCoef pfCounts_device[],	// the counts of all elements
	#endif	// MOD-BY-LEETEN 2013/07/03-END
	const unsigned int uNrOfDims, 
	const unsigned int uNrOfElements,
	unsigned int puKeys_device[],		// output: the keys of all elements. The keys are composed of bin and local_subscripts
	#if	0	// MOD-BY-LEETEN 2013/07/03-FROM:
	#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
	float pfCoefs_device[]				// output: the wavelet projection of the current wavelet 
	// ADD-BY-LEETEN 03/29/2013-BEGIN
	#else	// #if	!WITH_DOUBLE_COEF
	double pfCoefs_device[]				// output: the wavelet projection of the current wavelet 
	#endif	// #if	!WITH_DOUBLE_COEF
	// ADD-BY-LEETEN 03/29/2013-END
	#else	// MOD-BY-LEETEN 2013/07/03-TO:
	CudaDWT::typeCoef pfCoefs_device[]				// output: the wavelet projection of the current wavelet 
	#endif	// MOD-BY-LEETEN 2013/07/03-END
)
{
	/*
	int iX = blockIdx.x * blockDim.x + threadIdx.x;
	int iYZ = blockIdx.y * blockDim.y + threadIdx.y;
	int iY = (iYZ % iVolumeHeight_const);
	int iZ = (iYZ / iVolumeHeight_const);
	*/
	unsigned int uBlockId = blockIdx.x + blockIdx.y * gridDim.x;
	unsigned int uElement = uBlockId * blockDim.x + threadIdx.x;

	if( uElement < uNrOfElements )
	{
		uint4 u4BinSub = pu4BinSub_device[uElement];
		#if	0	// MOD-BY-LEETEN 2013/07/03-FROM:
		#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
		float fCount = pfCounts_device[uElement];
		// ADD-BY-LEETEN 03/29/2013-BEGIN
		#else	// #if	!WITH_DOUBLE_COEF
		double fCount = pfCounts_device[uElement];
		#endif	// #if	!WITH_DOUBLE_COEF	
		// ADD-BY-LEETEN 03/29/2013-END
		#else	// MOD-BY-LEETEN 2013/07/03-TO:
		CudaDWT::typeCoef fCount = pfCounts_device[uElement];
		#endif	// MOD-BY-LEETEN 2013/07/03-END

		unsigned int uKey = u4BinSub.x;
		unsigned int* puSub = &u4BinSub.y;
		int iWavelet = 1;
		for(unsigned int d = 0; d < uNrOfDims; d++)
		{
			unsigned int uSub = puSub[d];
			unsigned int w = puWaveletLengths_const[d];
			unsigned int l = puLevels_const[d];
			unsigned int uSubInWavelet = uSub % w;
			int iLocalWavelet = 1;
			if( 0 == l )
				iLocalWavelet = (int)w / 2 - (int)uSubInWavelet;
			else
			{
				if( uSubInWavelet < w / 2)
					iLocalWavelet = (int)uSubInWavelet;
				else
					iLocalWavelet = (int)w - (int)uSubInWavelet;
				iLocalWavelet *= -1;		
			}
			iWavelet *= iLocalWavelet;

			uKey *= puCoefLengths_const[d]/2;
			uKey += uSub / w;
		}

		#if	0	// MOD-BY-LEETEN 2013/07/03-FROM:
		#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
		pfCoefs_device[uElement] = fCount * (float)iWavelet;
		// ADD-BY-LEETEN 03/29/2013-BEGIN
		#else	// #if	!WITH_DOUBLE_COEF
		pfCoefs_device[uElement] = fCount * (double)iWavelet;
		#endif	// #if	!WITH_DOUBLE_COEF
		// ADD-BY-LEETEN 03/29/2013-END
		#else	// MOD-BY-LEETEN 2013/07/03-TO:
		pfCoefs_device[uElement] = fCount * (CudaDWT::typeCoef)iWavelet;
		#endif	// MOD-BY-LEETEN 2013/07/03-END
		puKeys_device[uElement] = uKey;
	}
}

