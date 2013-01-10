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
	unsigned int uElement = blockIdx.x * blockDim.x + threadIdx.x;

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

