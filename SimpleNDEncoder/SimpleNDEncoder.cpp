#include <math.h>
#include <assert.h>
#include <stdlib.h> 
#include "libclock.h"
#include "libopt.h"

#include "NrrdIO.h"
#include "SimpleND.h"

// namespace po = boost::program_options;

double dValueMin = HUGE_VAL;
double dValueMax = -HUGE_VAL;
Nrrd *nin;
CSimpleND<double> cSimpleND;
vector<double> vdData;

//! Convert the volume to an array of double type
template<typename T>
void 
_ConvertVolume
(
	const Nrrd *nin,
	double& dValueMin,
	double& dValueMax,
	vector<double>& vdData
)
{	
	T *data = (T*)nin->data;	

	// search for the range
	dValueMin = HUGE_VAL;
	dValueMax = -HUGE_VAL;
	for(int v = 0,	z = 0; z < (int)nin->axis[2].size; z++)
		for(int		y = 0; y < (int)nin->axis[1].size; y++)
			for(int x = 0; x < (int)nin->axis[0].size; x++, v++)
			{
				double dValue = (double)data[v];
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
	bool bIsPrintingTiming = true;
	LIBCLOCK_INIT(bIsPrintingTiming, __FUNCTION__);
	LIBCLOCK_BEGIN(bIsPrintingTiming);

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

	bool bIsOptParsed = BOPTParse(argv, argn, 1);

	/*
	int iSizeOfFullArrays = 0;
	int iNetCDFDeflateLevel = 0;
	int iNrOfBins = 8;	// iValueMax;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("size-of-full-arrays",		po::value<int>(&iSizeOfFullArrays)->default_value(iSizeOfFullArrays), 
			"Size (in MB) of the full arrays from all bin SATs")
		("netcdf-deflate-level",	po::value<int>(&iNetCDFDeflateLevel)->default_value(iNetCDFDeflateLevel), 
			"Deflate level for NetCDF file. The value is between 0 (store only) and 9 (maximal).")
		("n-bins",			po::value<int>(&iNrOfBins)->default_value(iNrOfBins), 
			"#Bins in the integral histograms")
		("vol-filepath",		po::value< string >(), 
			"The path to the input volum")
		("nc-filepath-prefix",		po::value< string >(), 
			"The prefixe of the output file");

	po::variables_map vm;
	po::store(po::parse_command_line(argn, argv, desc), vm);
	po::notify(vm);    

	if (vm.count("help")) {
		cout << desc << "\n";
		return 1;
	}
	const char *szVolFilePath = vm["vol-filepath"].as< string >().c_str();
	const char* szNcFilePathPrefix = vm["nc-filepath-prefix"].as< string >().c_str();
	*/
    
	assert(bIsOptParsed);
	assert(szVolFilePath);
	assert(szNcFilePathPrefix);

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

	// Step 1: Setup up the data size
	vector<size_t> vuDimLengths;
	size_t uNrOfValues = 1;	// ADD-BY-LEETEN 09/07/2012
	for(size_t d = 0; d < uNrOfDims; d++)
	{
		size_t uDimLength = (size_t)nin->axis[d].size;
		vuDimLengths.push_back( uDimLength );
		uNrOfValues *= uDimLength;
	}
	// MOD-BY-LEETEN 01/03/2013-FROM:	cSimpleND._Set(vuDimLengths, uNrOfBins);
	cSimpleND._Set(vuDimLengths, (WaveletSAT::typeBin)uNrOfBins);
	// MOD-BY-LEETEN 01/03/2013-END

	LOG_VAR(uNrOfBins);
	for(size_t d = 0; d < uNrOfDims; d++)
		LOG_VAR(vuDimLengths[d]);

	// Step 2: Allocate the needed #SATs
	// MOD-BY-LEETEN 01/03/2013-FROM:	cSimpleND._SetHistogram(uNrOfBins, dValueMin, dValueMax);
	#if 0 // MOD-BY-LEETEN 01/04/2013-FROM:
	cSimpleND._SetHistogram((WaveletSAT::typeBin)uNrOfBins, dValueMin, dValueMax);
	#else // MOD-BY-LEETEN 01/04/2013-TO:
	cSimpleND._SetHistogram(dValueMin, dValueMax);
	#endif // MOD-BY-LEETEN 01/04/2013-END
	// MOD-BY-LEETEN 01/03/2013-END
	cSimpleND._Allocate();
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
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	cSimpleND._Finalize();
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_BEGIN(bIsPrintingTiming);
	cSimpleND._SaveFile(szNcFilePathPrefix);
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_PRINT(bIsPrintingTiming);
	return 0;
}
