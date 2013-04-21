#include <math.h>
#include <assert.h>
#include <stdlib.h> 
#include "libclock.h"
#include "libopt.h"

#include "NrrdIO.h"
#include "SimpleND.h"

// namespace po = boost::program_options;

Nrrd *nin;
#if	!WITH_DOUBLE_COEF	// ADD-BY-LEETEN 03/29/2013
typedef float typeData;
// ADD-BY-LEETEN 03/29/2013-BEGIN
#else	// #if	!WITH_DOUBLE_COEF
typedef double typeData;
#endif	// #if	!WITH_DOUBLE_COEF
// ADD-BY-LEETEN 03/29/2013-END
typeData	dValueMin = (typeData)HUGE_VAL;
typeData	dValueMax = (typeData)-HUGE_VAL;
CSimpleND<typeData, typeData, WaveletSAT::typeBin, typeData> cSimpleND;
vector<typeData> vdData;

//! Convert the volume to an array of double type
template<typename T>
void 
_ConvertVolume
(
	const Nrrd *nin,
	typeData& dValueMin,
	typeData& dValueMax,
	vector<typeData>& vdData
)
{	
	T *data = (T*)nin->data;	

	// search for the range
	dValueMin = (typeData)HUGE_VAL;
	dValueMax = (typeData)-HUGE_VAL;
	for(int v = 0,	z = 0; z < (int)nin->axis[2].size; z++)
		for(int		y = 0; y < (int)nin->axis[1].size; y++)
			for(int x = 0; x < (int)nin->axis[0].size; x++, v++)
			{
				typeData dValue = (typeData)data[v];
				vdData.push_back(dValue);
				dValueMin = min(dValueMin, dValue);
				dValueMax = max(dValueMax, dValue);
			}
	LOG_VAR(dValueMin);
	LOG_VAR(dValueMax);
}

void
_ReadVolume
(
	const char* szPathFilename
)
{
	/* create a nrrd; at this point this is just an empty container */
	nin = nrrdNew();

	/* tell nrrdLoad to only read the header, not the data */
	NrrdIoState *nio = nrrdIoStateNew();
	nrrdIoStateSet(nio, nrrdIoStateSkipData, AIR_TRUE);

	/* read in the nrrd from file */
	if (nrrdLoad(nin, szPathFilename, nio)) {
		char *err = biffGetDone(NRRD);
		LOG_ERROR(fprintf(stderr, "%s", err));
		free(err);
		return;
	}

	/* we're done with the nrrdIoState, this sets it to NULL */
	nio = nrrdIoStateNix(nio);

	LOG_VAR(nrrdElementNumber(nin));
	LOG_VAR(nrrdElementSize(nin));
	nin->data = calloc(nrrdElementNumber(nin), nrrdElementSize(nin));

	if (nrrdLoad(nin, szPathFilename, NULL)) {
		char *err = biffGetDone(NRRD);
		LOG_ERROR(fprintf(stderr, "%s", err));
		free(err);
		return;
	}

	switch(nin->type)
	{
	case nrrdTypeChar:	_ConvertVolume<char>(nin, dValueMin, dValueMax, vdData);		break;
	case nrrdTypeUChar:	_ConvertVolume<unsigned char>(nin, dValueMin, dValueMax, vdData);	break;
	case nrrdTypeShort:	_ConvertVolume<short>(nin, dValueMin, dValueMax, vdData);		break;
	case nrrdTypeUShort:	_ConvertVolume<unsigned short>(nin, dValueMin, dValueMax, vdData);	break;
	case nrrdTypeInt:	_ConvertVolume<int>(nin, dValueMin, dValueMax, vdData);			break;
	case nrrdTypeUInt:	_ConvertVolume<unsigned int>(nin, dValueMin, dValueMax, vdData);	break;
	case nrrdTypeFloat:	_ConvertVolume<float>(nin, dValueMin, dValueMax, vdData);		break;

	default:
		break;
	}
}

int
main(int argn, char* argv[])
{
	_OPTInit();			// initialize the option parser

	char *szVolFilePath = NULL;
	_OPTAddStringVector(
		"--vol-filepath", 1,
		&szVolFilePath, szVolFilePath);

	int iNrOfBins = 8;	// iValueMax;
	_OPTAddIntegerVector(
		"--n-bins", 1,
		&iNrOfBins, iNrOfBins);

	char* szNcFilePathPrefix = NULL;
	_OPTAddStringVector(
		"--nc-filepath-prefix", 1,
		&szNcFilePathPrefix, szNcFilePathPrefix);

	int iIsVerbose = 0; 
	_OPTAddBoolean(
		"--is-verbose", &iIsVerbose, iIsVerbose);

	int iSizeOfFullArrays = 0;
	_OPTAddIntegerVector(
		"--size-of-full-arrays", 1,
		&iSizeOfFullArrays, iSizeOfFullArrays);
	_OPTAddComment("--size-of-full-arrays", 
		"Size (in MB) of the full arrays from all bin SATs");

	int iNetCDFDeflateLevel = 0;
	_OPTAddIntegerVector(
		"--netcdf-deflate-level", 1,
		&iNetCDFDeflateLevel, iNetCDFDeflateLevel);
	_OPTAddComment("--netcdf-deflate-level",
		"Deflate level for NetCDF file. The value is between 0 (store only) and 9 (maximal).");

	// ADD-BY-LEETEN 01/10/2013-BEGIN
	int iIsUsingGPUs = 0;
	_OPTAddBoolean("--is-using-gpus", &iIsUsingGPUs, iIsUsingGPUs);
	_OPTAddComment("--is-using-gpus",
		"The flag whether GPUs are used.");
	// ADD-BY-LEETEN 01/10/2013-END
	// ADD-BY-LEETEN 01/11/2013-BEGIN
	int iMaxNrOfEntriesOnGPUs = 4096;
	_OPTAddIntegerVector(
		"--max-n-entries-on-gpus", 1,
		&iMaxNrOfEntriesOnGPUs, iMaxNrOfEntriesOnGPUs);
	_OPTAddComment("--max-n-entries-on-gpus",
		"Max #Entries to be executed per GPU call. The unit is 1024.");
	// ADD-BY-LEETEN 01/11/2013-END
	// ADD-BY-LEETEN 01/11/2013-BEGIN
	int iTimingPrintingLevel = 1;
	_OPTAddIntegerVector(
		"--timing-printing-level", 1,
		&iTimingPrintingLevel, iTimingPrintingLevel);
	_OPTAddComment("--timing-printing-level",
		"The level to print the performance timing.");
	// ADD-BY-LEETEN 01/11/2013-END

	// ADD-BY-LEETEN 03/28/2013-BEGIN
	int iIsCompBinsOnly = 0;
	_OPTAddBoolean("--is-comp-bins-only", &iIsCompBinsOnly, iIsCompBinsOnly);
	// ADD-BY-LEETEN 03/28/2013-END

	// ADD-BY-LEETEN 04/20/2013-BEGIN
	int iIsUsingContourSpectrum = 0;
	_OPTAddBoolean(
		"--is-using-contour-spectrum", &iIsUsingContourSpectrum, iIsUsingContourSpectrum);
	_OPTAddComment("--is-using-contour-spectrum", 
		"Is the contour spectrum algorithm enabled?");
	// ADD-BY-LEETEN 04/20/2013-END

	bool bIsOptParsed = BOPTParse(argv, argn, 1);

	assert(bIsOptParsed);
	assert(szVolFilePath);
	assert(szNcFilePathPrefix);

	// ADD-BY-LEETEN 01/11/2013-BEGIN
	bool bIsPrintingTiming = (iTimingPrintingLevel>0)?true:false;
	LIBCLOCK_INIT(bIsPrintingTiming, __FUNCTION__);
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	// ADD-BY-LEETEN 01/11/2013-END

	LOG_VAR(szVolFilePath);
	_ReadVolume(szVolFilePath);

	size_t uNrOfDims = (size_t)nin->dim;
	size_t uMaxLevel = 0;
	size_t uWinSize = 1;
	size_t uNrOfBins = (size_t)iNrOfBins;

	LOG_VAR(iSizeOfFullArrays);	// ADD-BY-LEETEN 11/14/2012

	cSimpleND._SetInteger(CSimpleND<double>::SIZE_OF_FULL_ARRAYS, (long)iSizeOfFullArrays);
	#if WITH_NETCDF // ADD-BY-LEETEN 11/09/2012
	cSimpleND._SetInteger(CSimpleND<double>::DEFLATE_LEVEL, (long)iNetCDFDeflateLevel);
        #endif // #if WITH_NETCDF // ADD-BY-LEETEN 11/09/2012

	// ADD-BY-LEETEN 01/10/2012-BEGIN
	#if	WITH_CUDA
	cSimpleND._SetInteger(cSimpleND.IS_USING_GPUS, (long)iIsUsingGPUs);
	cSimpleND._SetInteger(cSimpleND.TIMING_PRINTING_LEVEL, (long)iTimingPrintingLevel - 1);	// ADD-BY-LEETEN 01/11/2013
	// ADD-BY-LEETEN 01/11/2013-BEGIN
	cSimpleND._SetInteger(cSimpleND.MAX_NR_OF_ELEMENTS_ON_THE_DEVICE, iMaxNrOfEntriesOnGPUs * 1024);
	// ADD-BY-LEETEN 01/11/2013-END
	#endif	// #if	WITH_CUDA
	// ADD-BY-LEETEN 01/10/2012-END

	// Step 1: Setup up the data size
	vector<size_t> vuDimLengths;
	size_t uNrOfValues = 1;	// ADD-BY-LEETEN 09/07/2012
	for(size_t d = 0; d < uNrOfDims; d++)
	{
		size_t uDimLength = (size_t)nin->axis[d].size;
		vuDimLengths.push_back( uDimLength );
		uNrOfValues *= uDimLength;
	}
	cSimpleND._Set(vuDimLengths, (WaveletSAT::typeBin)uNrOfBins);

	LOG_VAR(uNrOfBins);
	for(size_t d = 0; d < uNrOfDims; d++)
		LOG_VAR(vuDimLengths[d]);

	// Step 2: Allocate the needed #SATs
	cSimpleND._SetHistogram(dValueMin, dValueMax);
	cSimpleND._Allocate();
	// ADD-BY-LEETEN 04/20/2013-BEGIN
	cSimpleND._SetData(&vdData);
	cSimpleND._SetInteger(CSimpleND<double>::WITH_CONTOUR_SPECTRUM, iIsUsingContourSpectrum);
	// ADD-BY-LEETEN 04/20/2013-END
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_BEGIN(bIsPrintingTiming);
	// Step 3: Add the value to the SAT
	for(size_t i = 0; i < uNrOfValues; i++)
	{
		vector<size_t> vuPos;
		for(size_t 
			d = 0, uCoord = i; 
			d < nin->dim; 
			uCoord /= (size_t)nin->axis[d].size, d++)
			vuPos.push_back(uCoord % (size_t)nin->axis[d].size);

		cSimpleND._AddValue(vuPos, vdData[i]);
	}
	LIBCLOCK_END(bIsPrintingTiming);

	// Step 4: Finalize the SAT computation
	// ADD-BY-LEETEN 03/28/2013-BEGIN
	#if	WITH_SAT_FILE
	if( iIsCompBinsOnly )
	{
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		cSimpleND._SaveBins(szNcFilePathPrefix);
		LIBCLOCK_END(bIsPrintingTiming);
	}
	else
	{
	#endif	// #if	WITH_SAT_FILE
	// ADD-BY-LEETEN 03/28/2013-END
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	cSimpleND._Finalize();
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_BEGIN(bIsPrintingTiming);
	cSimpleND._SaveFile(szNcFilePathPrefix);
	LIBCLOCK_END(bIsPrintingTiming);

	// ADD-BY-LEETEN 03/28/2013-BEGIN
	#if	WITH_SAT_FILE
	}
	#endif	// #if	WITH_SAT_FILE
	// ADD-BY-LEETEN 03/28/2013-END

	LIBCLOCK_PRINT(bIsPrintingTiming);
	return 0;
}
