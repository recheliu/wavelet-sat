// ADD-BY-LEETEN 12/30/2012-BEGIN
#if		WITH_BOOST
#include <boost/filesystem.hpp>
using namespace boost::filesystem;
#endif	//	#if	WITH_BOOST
// ADD-BY-LEETEN 12/30/2012-END

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
// MOD-BY-LEETEN 12/29/2012-FROM:	CSimpleNDFile<double> cSimpleNDFile;
CSimpleNDFile<double, double> cSimpleNDFile;
// MOD-BY-LEETEN 12/29/2012-END
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

	// ADD-BY-LEETEN 12/30/2012-BEGIN
	int iEntropyWinRadius = 1;
	_OPTAddIntegerVector(
		"--entropy-win-radius", 1,
		&iEntropyWinRadius, iEntropyWinRadius);
	_OPTAddComment("--entropy-win-radius", 
		"Window Radius for entropy field computation");

	char* szEntropyFilepathPrefix = NULL;
	_OPTAddStringVector(
		"--entropy-filepath-prefix", 1,
		&szEntropyFilepathPrefix, szEntropyFilepathPrefix);
	_OPTAddComment("--entropy-filepath-prefix", 
		"Filepath prefix of the entropy field");
	// ADD-BY-LEETEN 12/30/2012-END

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
	LOG_VAR(szNcFilePath);	// ADD-BY-LEETEN 12/25/2012

	// load the WaveletSAT
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	LOG_VAR(iSizeOfFullArrays);	// ADD-BY-LEETEN 11/14/2012
	// MOD-BY-LEETEN 12/29/2012-FROM:	cSimpleNDFile._SetInteger(CSimpleNDFile<double>::SIZE_OF_FULL_ARRAYS, (long)iSizeOfFullArrays);
	cSimpleNDFile._SetInteger(cSimpleNDFile.SIZE_OF_FULL_ARRAYS, (long)iSizeOfFullArrays);
	// MOD-BY-LEETEN 12/29/2012-END
	cSimpleNDFile._LoadFile(szNcFilePath);
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_BEGIN(bIsPrintingTiming);
	cSimpleNDFile._ShowStatistics();
	LIBCLOCK_END(bIsPrintingTiming);

	#if	0	// DEL-BY-LEETEN 12/30/2012-BEGIN
	// ADD-By-LEETEN  12/29/2012-BEGIN
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	#if	0	// MOD-BY-LEETEN 12/30/2012-FROM:
	cSimpleNDFile._SetBoolean(cSimpleNDFile.PRINT_DECODE_BIN_TIMING, true);
	cSimpleNDFile._DecodeSAT();
	#else	// MOD-BY-LEETEN 12/30/2012-TO:
	cSimpleNDFile._SetInteger(cSimpleNDFile.PRINT_DECODE_BIN_TIMING, false);
	size_t uWinSize = 6;
	vector<int> viLeft, viRight;
	for(size_t d = 0; d < cSimpleNDFile.UGetNrOfDims(); d++)
	{
		viLeft.push_back(-(int)uWinSize);
		viRight.push_back(+(int)uWinSize);
	}
	cSimpleNDFile._ComputeEntropy(viLeft, viRight);
	#endif	// MOD-BY-LEETEN 12/30/2012-END
	LIBCLOCK_END(bIsPrintingTiming);	
	// ADD-By-LEETEN  12/29/2012-END
	#endif	// DEL-BY-LEETEN 12/30/2012-END

	if(iIsTestingQuery)
	{
		LIBCLOCK_BEGIN(bIsPrintingTiming);	// ADD-BY-LEETEN 12/30/2012

		// ADD-BY-LEETEN 12/28/2012-BEGIN
		// MOD-BY-LEETEN 12/28/2012-FROM:	cSimpleNDFile._SetBoolean(WaveletSAT::CSATSepDWTOutOfCore<double>::RESET_IO_COUNTERS, 0);
		// MOD-BY-LEETEN 12/30/2012-FROM:	cSimpleNDFile._SetBoolean(cSimpleNDFile.RESET_IO_COUNTERS, 0);
		cSimpleNDFile._SetInteger(cSimpleNDFile.RESET_IO_COUNTERS, 0);
		// MOD-BY-LEETEN 12/30/2012-END
		// MOD-BY-LEETEN 12/28/2012-END
		// ADD-BY-LEETEN 12/28/2012-END
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

		// DEL-BY-LEETEN 12/30/2012:	LIBCLOCK_BEGIN(bIsPrintingTiming);

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
		// ADD-BY-LEETEN 12/28/2012-BEGIN
		#if	0	// MOD-BY-LEETEN 12/28/2012-FROM:
		long lAccumNrOfIORequests;	cSimpleNDFile._GetInteger(WaveletSAT::CSATSepDWTOutOfCore<double>::ACCUM_NR_OF_IO_REQUESTS,	&lAccumNrOfIORequests);	LOG_VAR(lAccumNrOfIORequests);
		long lMaxNrOfIORequests;		cSimpleNDFile._GetInteger(WaveletSAT::CSATSepDWTOutOfCore<double>::MAX_NR_OF_IO_REQUESTS,		&lMaxNrOfIORequests);	LOG_VAR(lMaxNrOfIORequests);
		long lMinNrOfIORequests;		cSimpleNDFile._GetInteger(WaveletSAT::CSATSepDWTOutOfCore<double>::MIN_NR_OF_IO_REQUESTS,		&lMinNrOfIORequests);	LOG_VAR(lMinNrOfIORequests);
		#else	// MOD-BY-LEETEN 12/28/2012-TO:
		long lAccumNrOfIORequests;		cSimpleNDFile._GetInteger(cSimpleNDFile.ACCUM_NR_OF_IO_REQUESTS,	&lAccumNrOfIORequests);	LOG_VAR(lAccumNrOfIORequests);
		long lMaxNrOfIORequests;		cSimpleNDFile._GetInteger(cSimpleNDFile.MAX_NR_OF_IO_REQUESTS,		&lMaxNrOfIORequests);	LOG_VAR(lMaxNrOfIORequests);
		long lMinNrOfIORequests;		cSimpleNDFile._GetInteger(cSimpleNDFile.MIN_NR_OF_IO_REQUESTS,		&lMinNrOfIORequests);	LOG_VAR(lMinNrOfIORequests);
		#endif	// MOD-BY-LEETEN 12/28/2012-END
		// ADD-BY-LEETEN 12/28/2012-END

		// ADD-BY-LEETEN 12/30/2012-BEGIN
		/////////////////////////////////////////////////////////////
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		cSimpleNDFile._SetInteger(cSimpleNDFile.PRINT_DECODE_BIN_TIMING, false);
		vector<int> viLeft, viRight;
		for(size_t d = 0; d < cSimpleNDFile.UGetNrOfDims(); d++)
		{
			viLeft.push_back(-iEntropyWinRadius);
			viRight.push_back(+iEntropyWinRadius);
		}
		valarray<double> valEntropyField;
		cSimpleNDFile._ComputeEntropy(viLeft, viRight, valEntropyField);
		LIBCLOCK_END(bIsPrintingTiming);	

		/////////////////////////////////////////////////////////////
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		// Setup the file name
		ASSERT_OR_LOG(NULL != szEntropyFilepathPrefix, "");
		char szEntropyNhdrFilepath[NC_MAX_NAME+1];		sprintf(szEntropyNhdrFilepath,		"%s.nhdr",			szEntropyFilepathPrefix);
		char szEntropyRawFilepath[NC_MAX_NAME+1];		sprintf(szEntropyRawFilepath,		"%s.raw",			szEntropyFilepathPrefix);
		#if		WITH_BOOST
		path pathNhdr(szEntropyNhdrFilepath);
		path pathNhdrLeaf = pathNhdr.leaf();
		path pathNhdrDir =	pathNhdr.branch_path();
		path pathRaw(szEntropyRawFilepath);
		path pathRawLeaf = pathRaw.leaf();
		strcpy(szEntropyNhdrFilepath, pathNhdrLeaf.string().c_str());
		#else	// #if	WITH_BOOST
		#endif	// #if	WITH_BOOST

		//////////////////////////////////////////////
		// save the entropy field in Nrrd format
		#if	!WITH_SMART_PTR	// ADD-BY-LEETEN 12/30/2012
		TBuffer<float> pfEntropyField;		pfEntropyField.alloc(valEntropyField.size());
		// ADD-BY-LEETEN 12/30/2012-BEGIN
		#else	// #if	!WITH_SMART_PTR
		boost::shared_array<float> pfEntropyField(new float[valEntropyField.size()]);
		#endif	// #if	!WITH_SMART_PTR
		// ADD-BY-LEETEN 12/30/2012-END
		for(size_t d = 0; d < valEntropyField.size(); d++)
			pfEntropyField[d] = (float)valEntropyField[d];
		#if	!WITH_SMART_PTR		// ADD-BY-LEETEN 12/30/2012
		TBuffer<size_t> puSize;				puSize.alloc(cSimpleNDFile.UGetNrOfDims());
		// ADD-BY-LEETEN 12/30/2012-BEGIN
		#else	// #if	!WITH_SMART_PTR	
		boost::shared_array<size_t> puSize(new size_t[cSimpleNDFile.UGetNrOfDims()]);
		#endif	// #if	!WITH_SMART_PTR	
		// ADD-BY-LEETEN 12/30/2012-END
		for(size_t d = 0; d < cSimpleNDFile.UGetNrOfDims(); d++)
			puSize[d] = nin->axis[d].size;

		Nrrd *nrrdOut = nrrdNew();
		nrrdWrap_nva(nrrdOut, &pfEntropyField[0], nrrdTypeFloat, cSimpleNDFile.UGetNrOfDims(), &puSize[0]);
		nrrdSave(szEntropyNhdrFilepath, nrrdOut, NULL);

		// nrrdIoStateNix(nioOut);
		nrrdNix(nrrdOut);

		// now move the .nhdr and .raw to the destination folder
		#if		WITH_BOOST
		remove(pathNhdr);
		remove(pathRaw);
		rename(pathNhdrLeaf,	pathNhdr);
		rename(pathRawLeaf,		pathRaw);
		#else	// #if	WITH_BOOST
		#endif	// #if	WITH_BOOST
		LOG_VAR(szEntropyNhdrFilepath);	// ADD-BY-LEETEN 12/30/2012
		LIBCLOCK_END(bIsPrintingTiming);
		// ADD-BY-LEETEN 12/30/2012-END
	}

	LIBCLOCK_PRINT(bIsPrintingTiming);
	return 0;
}
