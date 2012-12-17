#include <math.h>
#include <assert.h>
#include <stdlib.h> 
#include "libclock.h"
#include "libopt.h"

#include "NrrdIO.h"
#include "SimpleNDFile.h"

double dValueMin = HUGE_VAL;
double dValueMax = -HUGE_VAL;
Nrrd *nin;
CSimpleNDFile<double> cSimpleNDFile;
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

	char* szNcFilePath = NULL;
	_OPTAddStringVector(
		"--nc-filepath", 1,
		&szNcFilePath, szNcFilePath);

	char *szVolFilePath = NULL;
	_OPTAddStringVector(
		"--vol-filepath", 1,
		&szVolFilePath, szVolFilePath);

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

	bool bIsOptParsed = BOPTParse(argv, argn, 1);
	assert(bIsOptParsed);
	assert(szVolFilePath);
	assert(szNcFilePath);

	// load the WaveletSAT
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	LOG_VAR(iSizeOfFullArrays);	// ADD-BY-LEETEN 11/14/2012
	cSimpleNDFile._SetLong(CSimpleNDFile<double>::SIZE_OF_FULL_ARRAYS, (long)iSizeOfFullArrays);
	cSimpleNDFile._LoadFile(szNcFilePath);
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_BEGIN(bIsPrintingTiming);
	cSimpleNDFile._ShowStatistics();
	LIBCLOCK_END(bIsPrintingTiming);

	if(iIsTestingQuery)
	{
		////////////////////////////////////////////////////////////////////////////
		// Now we can start to query SATs
		// load the data for testing
		LOG_VAR(szVolFilePath);
		_ReadVolume(szVolFilePath);

		size_t uNrOfDims = (size_t)nin->dim;
		size_t uNrOfTestingValues = (size_t)iNrOfTestingValues;
		size_t uWinSize = 1;
		size_t uNrOfBins = cSimpleNDFile.UGetNrOfBins();	// it will be decided later

		vector<size_t> vuDimLengths;
		size_t uNrOfValues = 1;	// ADD-BY-LEETEN 09/07/2012
		for(size_t d = 0; d < uNrOfDims; d++)
		{
			size_t uDimLength = (size_t)nin->axis[d].size;
			vuDimLengths.push_back( uDimLength );
			uNrOfValues *= uDimLength;
		}

		cSimpleNDFile._SetHistogram(uNrOfBins, dValueMin, dValueMax);

		// decide the threshld to filter numerical error
		double dThreshold = cSimpleNDFile.DGetThreshold();

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
			vector< pair<size_t, double> > vuBins;
			cSimpleNDFile._MapValueToBins(vuBase, vdData[uIndex], vuBins);

			size_t uValueBin = vuBins[0].first;
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
				cSimpleNDFile._GetAllSums(vuPos, vdIH);

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
