__constant__ unsigned int puWaveletLengths_const[CudaDWT::GPU_MAX_NR_OF_DIMS];
__constant__ unsigned int puLevels_const[CudaDWT::GPU_MAX_NR_OF_DIMS];
__constant__ unsigned int puCoefLengths_const[CudaDWT::GPU_MAX_NR_OF_DIMS];	

__global__ 
void 
_ProjToWavelet_kernel
(
	const uint4 pu4BinSub_device[],	// the tuples of <bin, data_subscripts> of all elements
	const CudaDWT::typeValue	pfCounts_device[],	// the counts of all elements
	const unsigned int uNrOfBins,		
	const unsigned int uNrOfDims, 
	const unsigned int uNrOfElements,
	CudaDWT::typeKey puKeys_device[],		// output: the keys of all elements. The keys are composed of bin and local_subscripts
	CudaDWT::typeCoef pfCoefs_device[]				// output: the wavelet projection of the current wavelet 
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
		CudaDWT::typeValue	fCount = pfCounts_device[uElement];

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
					iLocalWavelet =		(int)uSubInWavelet;
				else
					iLocalWavelet = (int)w - (int)uSubInWavelet;
				iLocalWavelet *= -1;		
			}
			iWavelet *= iLocalWavelet;
		}

		CudaDWT::typeKey uKey = 0;
		for(unsigned int di = 0; di < uNrOfDims; di++)
		{
			unsigned int d = uNrOfDims - 1 - di;

			unsigned int uSub = puSub[d];
			unsigned int w = puWaveletLengths_const[d];

			uKey *= puCoefLengths_const[d]/2;
			uKey += uSub / w;
		}

		uKey *= uNrOfBins;
		uKey += u4BinSub.x;

		pfCoefs_device[uElement] = (CudaDWT::typeCoef)(fCount * (CudaDWT::typeValue)iWavelet);
		puKeys_device[uElement] = uKey;
	}
}

