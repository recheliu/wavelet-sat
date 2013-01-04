#include <math.h>
#include <assert.h>
#include <stdlib.h> 
#include "libclock.h"
#include "libopt.h"

#include "NrrdIO.h"
#include "SimpleND.h"

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

	int iNrOfTestingValues = 0;
	_OPTAddIntegerVector(
		"--n-testing-values", 1,
		&iNrOfTestingValues, iNrOfTestingValues);
	

	int iIsTestingQuery = 0; 
	_OPTAddBoolean(
		"--is-testing-query", &iIsTestingQuery, iIsTestingQuery);

	int iIsVerbose = 0; 
	_OPTAddBoolean(
		"--is-verbose", &iIsVerbose, iIsVerbose);

	int iSizeOfFullArrays = 0;
	_OPTAddIntegerVector(
		"--size-of-full-arrays", 1,
		&iSizeOfFullArrays, iSizeOfFullArrays);
	_OPTAddComment("--size-of-full-arrays", 
		"Size (in MB) of the full arrays from all bin SATs");

	// ADD-BY-LEETEN 11/09/2012-BEGIN
	int iNetCDFDeflateLevel = 0;
	_OPTAddIntegerVector(
		"--netcdf-deflate-level", 1,
		&iNetCDFDeflateLevel, iNetCDFDeflateLevel);
	_OPTAddComment("--netcdf-deflate-level",
		"Deflate level for NetCDF file. The value is between 0 (store only) and 9 (maximal).");
	// ADD-BY-LEETEN 11/09/2012-END

	bool bIsOptParsed = BOPTParse(argv, argn, 1);
	assert(bIsOptParsed);
	assert(szVolFilePath);

	LOG_VAR(szVolFilePath);
	_ReadVolume(szVolFilePath);

	size_t uNrOfDims = (size_t)nin->dim;
	size_t uMaxLevel = 0;
	size_t uWinSize = 1;
	size_t uNrOfBins = (size_t)iNrOfBins;
	size_t uNrOfTestingValues = (size_t)iNrOfTestingValues;


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

	// ADD-BY-LEETEN 10/18/2012-BEGIN
	LOG_VAR(uNrOfBins);
	for(size_t d = 0; d < uNrOfDims; d++)
		LOG_VAR(vuDimLengths[d]);
	// ADD-BY-LEETEN 10/18/2012-END

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

	// ADD-BY-LEETEN 12/12/2012-BEGIN
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	cSimpleND._SaveFile(szVolFilePath);
	LIBCLOCK_END(bIsPrintingTiming);
	// ADD-BY-LEETEN 12/12/2012-END

	////////////////////////////////////////////////////////////////////////////
	// Now we can start to query SATs
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	cSimpleND._ShowStatistics();
	LIBCLOCK_END(bIsPrintingTiming);

	if(iIsTestingQuery)
	{
		// decide the threshld to filter numerical error
		double dThreshold = cSimpleND.DGetThreshold();

		LIBCLOCK_BEGIN(bIsPrintingTiming);

		size_t uNrOfIHs = 1 << uNrOfDims;
		vector< vector<size_t> > vvuOffsets;
		vvuOffsets.resize(uNrOfIHs);
		for(size_t i = 0; i < uNrOfIHs; i++)
		{
			vvuOffsets[i].resize(uNrOfDims);
			for(size_t 
				d = 0, j = i; 
				d < uNrOfDims; 
				d++, j /= 2)
				vvuOffsets[i][d] = uWinSize * (j % 2);
		}
		LIBCLOCK_END(bIsPrintingTiming);


		LIBCLOCK_BEGIN(bIsPrintingTiming);
		for(size_t t = 0; t < uNrOfTestingValues; t++)
		{
			vector<size_t> vuBase;
			size_t uIndex = 0;
			if( iIsVerbose )
				printf("B(");
			for(size_t 
				d = 0, uDimLengthProduct = 1; 
				d < uNrOfDims; 
				uDimLengthProduct *= vuDimLengths[d], d++)
			{
				size_t uPos = uWinSize + rand() % (vuDimLengths[d] - uWinSize);
				vuBase.push_back(uPos);
				uIndex += uPos * uDimLengthProduct;

				if( iIsVerbose )
					printf("%3d,", (int)uPos);
			}
			// MOD-BY-LEETEN 01/03/2013-FROM:			vector< pair<size_t, double> > vuBins;
			vector< pair<WaveletSAT::typeBin, double> > vuBins;
			// MOD-BY-LEETEN 01/03/2013-END
			cSimpleND._MapValueToBins(vuBase, vdData[uIndex], vuBins);

			// MOD-BY-LEETEN 01/03/2013-FROM:			size_t uValueBin = vuBins[0].first;
			WaveletSAT::typeBin uValueBin = vuBins[0].first;
			// MOD-BY-LEETEN 01/03/2013-END
			if( iIsVerbose )
				printf(")=\t%d,\n", (int)uValueBin);	// vdData[uIndex]);	// 

			vector<double> vdH;
			vdH.resize(uNrOfBins);
			for(size_t i = 0; i < uNrOfIHs; i++)
			{
				vector<size_t> vuPos;
				int iSign = 1;
				for(size_t d = 0; d < uNrOfDims; d++)
				{
					vuPos.push_back(vuBase[d] - vvuOffsets[i][d]);
					iSign *= (vvuOffsets[i][d])?(-1):(+1);
				}
				vector<double> vdIH;
				cSimpleND._GetAllSums(vuPos, vdIH);

				for(size_t b = 0; b < uNrOfBins; b++)
					vdH[b] += iSign * vdIH[b]; 
			}

			double dError = 0.0;
			for(size_t b = 0; b < uNrOfBins; b++)
				if(b == uValueBin)
					dError += pow(1.0 - vdH[b], 2.0);
				else
					dError += pow(vdH[b], 2.0);

			if( iIsVerbose )
			{
				printf("H:");
				for(size_t b = 0; b < uNrOfBins; b++)
					if( fabs(vdH[b]) > dThreshold )
					  printf( "\t\t%d:%+.2f\n", (unsigned int)b, vdH[b]);
				printf("E:%f\n", dError);
			}
		}

		LIBCLOCK_END(bIsPrintingTiming);
	}

	LIBCLOCK_PRINT(bIsPrintingTiming);
	return 0;
}
