__global__ 
void 
_MarkSegments_kernel
(
	unsigned int puKey_devices[],
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
	// MOD-BY-LEETEN 01/11/2013-FROM:	unsigned int uElement = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int uBlockId = blockIdx.x + blockIdx.y * gridDim.x;
	unsigned int uElement = uBlockId * blockDim.x + threadIdx.x;
	// MOD-BY-LEETEN 01/11/2013-END

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
