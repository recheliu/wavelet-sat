__global__ 
void 
_MarkSegments_kernel
(
	// MOD-BY-LEETEN 2013/07/31-FROM:	unsigned int puKey_devices[],
	CudaDWT::typeKey puKey_devices[],
	// MOD-BY-LEETEN 2013/07/31-END
	unsigned int uNrOfElements,
	unsigned int puiSegFlags_device[]
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
		unsigned int uiSeg = 0;
		if( 0 == uElement )
			uiSeg = 1;
		else
			uiSeg = (puKey_devices[uElement] != puKey_devices[uElement - 1])?1:0;

		puiSegFlags_device[uElement] = uiSeg;
	}
}

