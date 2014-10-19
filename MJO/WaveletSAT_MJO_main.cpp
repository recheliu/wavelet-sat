#include <math.h>
#include <assert.h>
#include <stdlib.h> 
#include "libclock.h"
#include "libopt.h"

#include "SimpleND.h"

// namespace po = boost::program_options;

#if	!WITH_DOUBLE_COEF	
	typedef float typeData;
#else	// #if	!WITH_DOUBLE_COEF
	typedef double typeData;
#endif	// #if	!WITH_DOUBLE_COEF

typeData	dValueMin = (typeData)HUGE_VAL;
typeData	dValueMax = (typeData)-HUGE_VAL;
#if	!WITH_DOUBLE_COEF	
	CSimpleND<typeData, float, WaveletSAT::typeBin, float> cSimpleND;
#else	// #if	!WITH_DOUBLE_COEF	
	CSimpleND<typeData, typeData, WaveletSAT::typeBin, typeData> cSimpleND;
#endif	// #if	!WITH_DOUBLE_COEF

vector<typeData> vdData;
vector<size_t> vuDimLengths;
size_t uNrOfValues;

void
_ReadVolume
(
	const char *szVecDirPath,
	const char *szVecFileName,
	const char *szVecVarName,
	void* _Reserved = NULL
)
{
	char szVecFilePath[NC_MAX_NAME];
	sprintf(szVecFilePath, "%s/%s", szVecDirPath, szVecFileName);

	int ncId;
	ASSERT_NETCDF(nc_open(
    	szVecFilePath,
    	NC_NOWRITE,
    	&ncId));

	int var;
	ASSERT_NETCDF(nc_inq_varid(
		ncId, 
		szVecVarName, 
		&var) );

	int iNrOfDims;
	ASSERT_NETCDF(nc_inq_varndims(
		ncId, 
		var,
		&iNrOfDims));

	int pdimIDs[NC_MAX_DIMS];
	ASSERT_NETCDF(nc_inq_vardimid(
		ncId, 
		var,
		pdimIDs));

	vuDimLengths.clear();
	for(size_t d = 0; d < (size_t) iNrOfDims; d++)
	{
		size_t uDimLen;
		ASSERT_NETCDF(nc_inq_dimlen(
			ncId, 
			pdimIDs[iNrOfDims - 1 - d], 
			&uDimLen));
		vuDimLengths.push_back(uDimLen);
	}

	// set up the dim length for this slice
	size_t puStarts[NC_MAX_DIMS];
	size_t puCounts[NC_MAX_DIMS];
	for(size_t d = 0; d < (size_t)iNrOfDims; d++)
	{
		puStarts[d] = 0;
		puCounts[d] = vuDimLengths[iNrOfDims - 1 - d];
	}

	// only consider t = 0;
	vuDimLengths[iNrOfDims - 1] = puCounts[0] = 1;

	// only consider z = 12;
	vuDimLengths[iNrOfDims - 1 - 1] = puCounts[1] = 1;
	puStarts[1] = 12;


	uNrOfValues = 1;
	for(size_t d = 0; d < (size_t)iNrOfDims; d++)
		uNrOfValues *= (size_t)puCounts[d];

	vdData.resize(uNrOfValues);

	// load the data
	ASSERT_NETCDF(nc_get_vara(
		ncId, 
		var, 
		&puStarts[0],
		&puCounts[0], 
		vdData.data()));

	ASSERT_NETCDF(nc_close(ncId));

	for(size_t v = 0; v < uNrOfValues; v++)
	{
		dValueMin = min(dValueMin, vdData[v]);
		dValueMax = max(dValueMax, vdData[v]);
	}
	dValueMin = max(dValueMin, (typeData)0.0);

	// now remove the first two dimensions, which are z and t.
	vuDimLengths.pop_back();
	vuDimLengths.pop_back();
}

int
main(int argn, char* argv[])
{
	_OPTInit();			// initialize the option parser

	char *szDirPath = NULL;
	char *szFileName = NULL;
	char *szVarName = NULL;
	_OPTAddStringVector(
		"--data", 3,
			&szDirPath,		szDirPath,
			&szFileName,	szFileName,
			&szVarName,		szVarName
		);

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

	int iIsUsingContourSpectrum = 0;
	_OPTAddBoolean(
		"--is-using-contour-spectrum", &iIsUsingContourSpectrum, iIsUsingContourSpectrum);
	_OPTAddComment("--is-using-contour-spectrum", 
		"Is the contour spectrum algorithm enabled?");

	int iIsCompBinsOnly = 0;
	_OPTAddBoolean("--is-comp-bins-only", &iIsCompBinsOnly, iIsCompBinsOnly);

	bool bIsOptParsed = BOPTParse(argv, argn, 1);

	assert(bIsOptParsed);
	assert(szDirPath);
	assert(szFileName);
	assert(szVarName);
	assert(szNcFilePathPrefix);

	bool bIsPrintingTiming = (iTimingPrintingLevel>0)?true:false;
	LIBCLOCK_INIT(bIsPrintingTiming, __FUNCTION__);
	LIBCLOCK_BEGIN(bIsPrintingTiming);

	_ReadVolume(szDirPath, szFileName, szVarName);

	size_t uNrOfDims = vuDimLengths.size();
	size_t uMaxLevel = 0;
	size_t uWinSize = 1;
	size_t uNrOfBins = (size_t)iNrOfBins;

	LOG_VAR(iSizeOfFullArrays);	

	cSimpleND._SetInteger(CSimpleND<double>::SIZE_OF_FULL_ARRAYS, (long)iSizeOfFullArrays);
	#if WITH_NETCDF 
	cSimpleND._SetInteger(CSimpleND<double>::DEFLATE_LEVEL, (long)iNetCDFDeflateLevel);
    #endif // #if WITH_NETCDF 

	#if	WITH_CUDA
	cSimpleND._SetInteger(cSimpleND.IS_USING_GPUS, (long)iIsUsingGPUs);
	cSimpleND._SetInteger(cSimpleND.TIMING_PRINTING_LEVEL, (long)iTimingPrintingLevel - 1);	
	cSimpleND._SetInteger(cSimpleND.MAX_NR_OF_ELEMENTS_ON_THE_DEVICE, iMaxNrOfEntriesOnGPUs * 1024);
	#endif	// #if	WITH_CUDA

	// Step 1: Setup up the data size
	cSimpleND._Set(vuDimLengths, (WaveletSAT::typeBin)uNrOfBins);

	LOG_VAR(uNrOfBins);
	for(size_t d = 0; d < uNrOfDims; d++)
		LOG_VAR(vuDimLengths[d]);

	// Step 2: Allocate the needed #SATs
	cSimpleND._SetHistogram(dValueMin, dValueMax);
	cSimpleND._Allocate();
	cSimpleND._SetData(&vdData);
	cSimpleND._SetInteger(CSimpleND<double>::WITH_CONTOUR_SPECTRUM, iIsUsingContourSpectrum);
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_BEGIN(bIsPrintingTiming);
	// Step 3: Add the value to the SAT
	for(size_t i = 0; i < uNrOfValues; i++)
	{
		vector<size_t> vuPos;
		for(size_t 
			d = 0, uCoord = i; 
			d < vuDimLengths.size(); 
			uCoord /= (size_t)vuDimLengths[d], d++)
			vuPos.push_back(uCoord % (size_t)vuDimLengths[d]);

		cSimpleND._AddValue(vuPos, vdData[i]);
	}
	LIBCLOCK_END(bIsPrintingTiming);

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
		// Step 4: Finalize the SAT computation
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		cSimpleND._Finalize();
		LIBCLOCK_END(bIsPrintingTiming);

		LIBCLOCK_BEGIN(bIsPrintingTiming);
		cSimpleND._SaveFile(szNcFilePathPrefix);
		LIBCLOCK_END(bIsPrintingTiming);

	#if	WITH_SAT_FILE
	}
	#endif	// #if	WITH_SAT_FILE

	LIBCLOCK_PRINT(bIsPrintingTiming);
	return 0;
}
