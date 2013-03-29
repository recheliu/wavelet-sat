#include <math.h>
#include <assert.h>
#include <stdlib.h> 
#include "libclock.h"
#include "libopt.h"

#include "Vector2D.h"

CVector2D<> cVector2D;

int
main(int argn, char* argv[])
{
	_OPTInit();			// initialize the option parser

	char *szVecDirPath = NULL;
	char *szUVecFileName = NULL;
	char *szUVecVarName = NULL;
	char *szVVecFileName = NULL;
	char *szVVecVarName = NULL;
	_OPTAddStringVector(
		"--vec-filepaths", 5,
			&szVecDirPath, szVecDirPath,
			&szUVecFileName, szUVecFileName,
			&szUVecVarName, szUVecVarName,
			&szVVecFileName, szVVecFileName, 
			&szVVecVarName, szVVecVarName
		);
	_OPTAddNames(
		"--vec-filepaths", 5, 
			"VecDirPath",
			"UVecFileName", 
			"VVecFileName", 
			"UVecVarName", 
			"VVecVarName"
		);
	_OPTAddComment(
		"--vec-filepaths", 
			"The NetCDF files that store the U and V components of the vector field. The variable name is also specified"
		);

	int iDepth = 0;	// iValueMax;
	_OPTAddIntegerVector(
		"--depth", 1,
		&iDepth, iDepth);

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

	int iIsUsingGPUs = 0;
	_OPTAddBoolean("--is-using-gpus", &iIsUsingGPUs, iIsUsingGPUs);
	_OPTAddComment("--is-using-gpus",
		"The flag whether GPUs are used.");

	int iMaxNrOfEntriesOnGPUs = 4096;
	_OPTAddIntegerVector(
		"--max-n-entries-on-gpus", 1,
		&iMaxNrOfEntriesOnGPUs, iMaxNrOfEntriesOnGPUs);
	_OPTAddComment("--max-n-entries-on-gpus",
		"Max #Entries to be executed per GPU call. The unit is 1024.");

	int iTimingPrintingLevel = 1;
	_OPTAddIntegerVector(
		"--timing-printing-level", 1,
		&iTimingPrintingLevel, iTimingPrintingLevel);
	_OPTAddComment("--timing-printing-level",
		"The level to print the performance timing.");

	// ADD-BY-LEETEN 03/28/2013-BEGIN
	int iIsCompBinsOnly = 1;
	_OPTAddBoolean("--is-comp-bins-only", &iIsCompBinsOnly, iIsCompBinsOnly);
	// ADD-BY-LEETEN 03/28/2013-END

	bool bIsOptParsed = BOPTParse(argv, argn, 1);

	assert(bIsOptParsed);
	assert(szNcFilePathPrefix);

	bool bIsPrintingTiming = (iTimingPrintingLevel>0)?true:false;
	LIBCLOCK_INIT(bIsPrintingTiming, __FUNCTION__);
	LIBCLOCK_BEGIN(bIsPrintingTiming);

	cVector2D._SetInteger(cVector2D.SIZE_OF_FULL_ARRAYS,	(long)iSizeOfFullArrays);
	cVector2D._SetInteger(cVector2D.DEFLATE_LEVEL,			(long)iNetCDFDeflateLevel);
	cVector2D._SetInteger(cVector2D.NR_OF_BINS,				(long)iNrOfBins);
	cVector2D._SetInteger(cVector2D.DEPTH,					(long)iDepth);

	#if	WITH_CUDA
	cVector2D._SetInteger(cVector2D.IS_USING_GPUS, (long)iIsUsingGPUs);
	cVector2D._SetInteger(cVector2D.TIMING_PRINTING_LEVEL, (long)iTimingPrintingLevel - 1);	// ADD-BY-LEETEN 01/11/2013
	cVector2D._SetInteger(cVector2D.MAX_NR_OF_ELEMENTS_ON_THE_DEVICE, iMaxNrOfEntriesOnGPUs * 1024);
	#endif	// #if	WITH_CUDA

	#if	0	// DEL-BY-LEETEN 03/19/2013-BEGIN
	// ADD-By-LEETEN 02/19/2013-BEGIN
	LIBCLOCK_END(bIsPrintingTiming);
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	// ADD-By-LEETEN 02/19/2013-END
	#endif	// DEL-BY-LEETEN 03/19/2013-END
	cVector2D._Load(
		szVecDirPath,
		szUVecFileName,
		szUVecVarName,
		szVVecFileName,
		szVVecVarName);

	LIBCLOCK_END(bIsPrintingTiming);

	#if	0	// DEL-BY-LEETEN 03/28/2013-BEGIN
	// ADD-BY-LEETEN 03/19/2013-BEGIN
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	LIBCLOCK_END(bIsPrintingTiming);
	// ADD-BY-LEETEN 03/19/2013-END
	#endif	// DEL-BY-LEETEN 03/28/2013-END

	LIBCLOCK_BEGIN(bIsPrintingTiming);
	cVector2D._Encode();
	LIBCLOCK_END(bIsPrintingTiming);

	// ADD-BY-LEETEN 03/28/2013-BEGIN
	#if	WITH_SAT_FILE
	if( iIsCompBinsOnly )
	{
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		cVector2D._SaveBins(szNcFilePathPrefix);
		LIBCLOCK_END(bIsPrintingTiming);
	}
	else
	{
	#endif	// #if	WITH_SAT_FILE
	// ADD-BY-LEETEN 03/28/2013-END
	cVector2D._Finalize();	// ADD-BY-LEETEN 03/28/2013

	LIBCLOCK_BEGIN(bIsPrintingTiming);
	cVector2D._SaveFile(szNcFilePathPrefix);
	LIBCLOCK_END(bIsPrintingTiming);

	// ADD-BY-LEETEN 03/28/2013-BEGIN
	#if	WITH_SAT_FILE
	}
	#endif	// #if	WITH_SAT_FILE
	// ADD-BY-LEETEN 03/28/2013-END

	LIBCLOCK_PRINT(bIsPrintingTiming);
	return 0;
}
