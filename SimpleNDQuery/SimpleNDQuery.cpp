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
CSimpleNDFile<double, double, WaveletSAT::typeBin, double> cSimpleNDFile;
vector<double> vdData;

// ADD-BY-LEETEN 03/28/2013-BEGIN
#include "lognc.h"

vector<WaveletSAT::typeBin> vusCoefBins;
vector<WaveletSAT::typeSum> vdCoefSums;
vector<WaveletSAT::CSATFileNetCDF::typeCoefCount>		vusCoefCounts;
vector<WaveletSAT::CSATFileNetCDF::typeCoefOffset>		vuCoefOffsets;

template<typename T>
void
_ReadVar(
	int iNcId,
	const char *szVar,
	vector<T>& vData 
)
{
	int piDimIds[NC_MAX_DIMS];
	int iNrOfDims;
	int ncVarId;
	ASSERT_NETCDF(nc_inq_varid(iNcId, szVar, &ncVarId));
	ASSERT_NETCDF(nc_inq_varndims(iNcId, ncVarId, &iNrOfDims));
	ASSERT_NETCDF(nc_inq_vardimid(iNcId, ncVarId, piDimIds));
	size_t uDataLen = 1;
	for(size_t d = 0; d < (size_t)iNrOfDims; d++)
	{
		size_t uDimLen;
		ASSERT_NETCDF(nc_inq_dimlen(iNcId, piDimIds[d], &uDimLen));
		uDataLen *= uDimLen;
	}
	vData.resize(uDataLen);
	ASSERT_NETCDF(
		nc_get_var(iNcId, ncVarId, (void*)&vData.data()[0]));
}

void
_ReadBins(
	const char *szFilepath
)
{
	int iNcId;

	#if !WITH_NETCDF4 
	ASSERT_NETCDF(nc_open(
    		szFilepath,
    		NC_NOWRITE,
    		&iNcId));
	#else	// #if !WITH_NETCDF4
	ASSERT_NETCDF(nc_open(
    		szFilepath,
    		NC_NOWRITE | NC_NETCDF4,
    		&iNcId));
	#endif // #if !WITH_NETCDF4

	_ReadVar<WaveletSAT::typeBin>(iNcId, WaveletSAT::CSATFileNetCDF::SZGetVarCoefBin(), vusCoefBins);
	_ReadVar<WaveletSAT::typeSum>(iNcId, WaveletSAT::CSATFileNetCDF::SZGetVarCoefSum(), vdCoefSums);
	_ReadVar<WaveletSAT::CSATFileNetCDF::typeCoefCount>(iNcId,	WaveletSAT::CSATFileNetCDF::SZGetVarCoefCount(), vusCoefCounts);
	_ReadVar<WaveletSAT::CSATFileNetCDF::typeCoefOffset>(iNcId, WaveletSAT::CSATFileNetCDF::SZGetVarCoefOffset(), vuCoefOffsets);
	ASSERT_NETCDF(
		nc_close(iNcId) );
}
// ADD-BY-LEETEN 03/28/2013-END

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
	
	// ADD-BY-LEETEN 01/09/2013-BEGIN
	int iQueryWinLength = 1;
	_OPTAddIntegerVector(
		"--query-win-length", 1,
		&iQueryWinLength, iQueryWinLength);
	// ADD-BY-LEETEN 01/09/2013-END

	int iIsTestingQuery = 0; 
	_OPTAddBoolean(
		"--is-testing-query", &iIsTestingQuery, iIsTestingQuery);

	int iIsVerbose = 0; 
	_OPTAddBoolean(
		"--is-verbose", &iIsVerbose, iIsVerbose);

	// ADD-BY-LEETEN 12/30/2012-BEGIN
	// ADD-BY-LEETEN 12/31/2012-BEGIN
	int iIsComputingEntropy = 0; 
	_OPTAddBoolean(
		"--is-computing-entropy", &iIsComputingEntropy, iIsComputingEntropy);
	_OPTAddComment("--is-computing-entropy", 
		"The flag indicating whether the entropy field is computed or not.");
	// ADD-BY-LEETEN 12/31/2012-END

	int iEntropyWinRadius = 1;
	_OPTAddIntegerVector(
		"--entropy-win-radius", 1,
		&iEntropyWinRadius, iEntropyWinRadius);
	_OPTAddComment("--entropy-win-radius", 
		"Window Radius for entropy field computation");

	// ADD-BY-LEETEN 2013/12/01-BEGIN
	int iNrOfEntropyBins = 0;
	_OPTAddIntegerVector(
		"--n-entropy-bins", 1,
		&iNrOfEntropyBins, iNrOfEntropyBins);
	_OPTAddComment("--n-entropy-bins", 
		"#Bins for entropy field computation");
	// ADD-BY-LEETEN 2013/12/01-END

	char* szEntropyFilepathPrefix = NULL;
	_OPTAddStringVector(
		"--entropy-filepath-prefix", 1,
		&szEntropyFilepathPrefix, szEntropyFilepathPrefix);
	_OPTAddComment("--entropy-filepath-prefix", 
		"Filepath prefix of the entropy field");
	// ADD-BY-LEETEN 12/30/2012-END

	// ADD-BY-LEETEN 03/28/2013-BEGIN
	char *szBinFilePath = NULL;
	_OPTAddStringVector(
		"--bin-filepath", 1,
		&szBinFilePath, szBinFilePath);
	// ADD-BY-LEETEN 03/28/2013-END

	int iSizeOfFullArrays = 0;
	_OPTAddIntegerVector(
		"--size-of-full-arrays", 1,
		&iSizeOfFullArrays, iSizeOfFullArrays);
	_OPTAddComment("--size-of-full-arrays", 
		"Size (in MB) of the full arrays from all bin SATs");

	// ADD-BY-LEETEN 01/23/2013-BEGIN
	enum {
		STAT_MEAN,
		STAT_COUNT,
		STAT_STDDEV,
		STAT_ENTROPY,
		NR_OF_STATS,
		STAT_DEFAULT = STAT_MEAN 
	};
	int iStat = STAT_DEFAULT;
	_OPTAddEnum(
		"--stat", &iStat, iStat, NR_OF_STATS,
		"mean",		STAT_MEAN,
		"count",	STAT_COUNT,
		"stddev",	STAT_STDDEV,
		"entropy",	STAT_ENTROPY);

	int iBlockLevel = 2;
	_OPTAddIntegerVector(
		"--block-level", 1,
		&iBlockLevel, iBlockLevel);

	char *szStatFilepathPrefix = NULL;
	_OPTAddStringVector(
		"--stat-filepath-prefix", 1,
		&szStatFilepathPrefix, szStatFilepathPrefix);

	int iIsComputingBlockStat = 0; 
	_OPTAddBoolean(
		"--is-computing-block-stat", &iIsComputingBlockStat, iIsComputingBlockStat);
	// ADD-BY-LEETEN 01/23/2013-END

	// ADD-BY-LEETEN 2013/07/03-BEGIN
	int iIsTestingBruteForce = 0; 
	_OPTAddBoolean(
		"--is-testing-brute-force", &iIsTestingBruteForce, iIsTestingBruteForce);
	// ADD-BY-LEETEN 2013/07/03-END

	bool bIsOptParsed = BOPTParse(argv, argn, 1);
	assert(bIsOptParsed);
	assert(szNcFilePath);
	LOG_VAR(szNcFilePath);	// ADD-BY-LEETEN 12/25/2012

	// load the WaveletSAT
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	LOG_VAR(iSizeOfFullArrays);	// ADD-BY-LEETEN 11/14/2012
	cSimpleNDFile._SetInteger(cSimpleNDFile.SIZE_OF_FULL_ARRAYS, (long)iSizeOfFullArrays);
	cSimpleNDFile._LoadFile(szNcFilePath);
	LIBCLOCK_END(bIsPrintingTiming);

	if(iIsTestingQuery)
	{
		// ADD-BY-LEETEN 03/28/2013-BEGIN
		if( szBinFilePath )
		{
			LIBCLOCK_BEGIN(bIsPrintingTiming);	// ADD-BY-LEETEN 12/30/2012

			cSimpleNDFile._SetInteger(cSimpleNDFile.RESET_IO_COUNTERS, 0);

			////////////////////////////////////////////////////////////////////////////
			// Now we can start to query SATs
			// load the data for testing
			LOG_VAR(szBinFilePath);
			size_t uNrOfDims = cSimpleNDFile.UGetNrOfDims();
			size_t uNrOfTestingValues = (size_t)iNrOfTestingValues;
			size_t uWinSize = iQueryWinLength;
			size_t uNrOfBins = cSimpleNDFile.UGetNrOfBins();	// it will be decided later
			vector<size_t> vuDimLengths;
			cSimpleNDFile._GetDataSize(vuDimLengths);
			size_t uNrOfValues = 1;	// ADD-BY-LEETEN 09/07/2012
			for(size_t d = 0; d < uNrOfDims; d++)
				uNrOfValues *= vuDimLengths[d];

			// decide the threshld to filter numerical error
			double dThreshold = cSimpleNDFile.DGetThreshold();

			_ReadBins(szBinFilePath);

			LIBCLOCK_END(bIsPrintingTiming);

			LIBCLOCK_BEGIN(bIsPrintingTiming);

			vector<size_t> vuWinSize;
			size_t uWinLength = 1;
			for(size_t d = 0; d < uNrOfDims; d++)
			{
				uWinLength *= uWinSize;
				vuWinSize.push_back(uWinSize);
			}

			for(size_t t = 0; t < uNrOfTestingValues; t++)
			{
				vector<size_t> vuBase;
				size_t uIndex = 0;

				for(size_t 
					d = 0, uDimLengthProduct = 1; 
					d < uNrOfDims; 
					uDimLengthProduct *= vuDimLengths[d], d++)
				{
					if(uWinSize)
						vuBase.push_back(uWinSize + rand() % (vuDimLengths[d] - uWinSize));
					else
						vuBase.push_back(rand() % vuDimLengths[d]);
				}

				vector<WaveletSAT::typeSum> vdH;
				vdH.assign(uNrOfBins, (WaveletSAT::typeSum)0);

				vector<WaveletSAT::typeSum> vdRefH;
				vdRefH.assign(uNrOfBins, (WaveletSAT::typeSum)0);

				if(uWinSize)
				{
					vector<size_t> vuOffset;
					vuOffset.resize(uNrOfDims);
					for(size_t d = 0; d < uNrOfDims; d++)
						vuOffset[d] = vuBase[d] - uWinSize;
					// ADD-BY-LEETEN 2013/07/03-BEGIN
					if( !iIsTestingBruteForce )
					// ADD-BY-LEETEN 2013/07/03-END
					cSimpleNDFile._GetRegionSums(vuOffset, vuBase, vdH);

					for(size_t w = 0; w < uWinLength; w++)
					{
						vector<size_t> vuPosInWin;
						WaveletSAT::_ConvertIndexToSub(w, vuPosInWin, vuWinSize);

						vector<size_t> vuPos;
						vuPos.resize(uNrOfDims);
						for(size_t d = 0; d < uNrOfDims; d++)
							vuPos[d] = vuOffset[d] + 1 + vuPosInWin[d];

						size_t uIndex = WaveletSAT::UConvertSubToIndex(vuPos, vuDimLengths);
						size_t uOffset = vuCoefOffsets[uIndex];
						for(size_t bi = 0; bi < vusCoefCounts[uIndex]; bi++)
						{
							vdRefH[vusCoefBins[uOffset+bi]] += (WaveletSAT::typeSum)vdCoefSums[uOffset+bi];
						}
					}
				}
				else
				{
					cSimpleNDFile._GetAllSums(vuBase, vdH);

					vector<size_t> vuBaseLengths;
					vuBaseLengths.resize(uNrOfDims);
					size_t uBaseSize= 1;
					for(size_t d = 0; d < uNrOfDims; d++)
					{
						uBaseSize *= vuBase[d] + 1;
						vuBaseLengths[d] = vuBase[d] + 1;
					}

					vector<size_t> vuB;
					vuB.resize(uNrOfDims);
					for(size_t b = 0; b < uBaseSize; b++)
					{
						WaveletSAT::_ConvertIndexToSub(b, vuB, vuBaseLengths);
						size_t uIndex = WaveletSAT::UConvertSubToIndex(vuB, vuDimLengths);
						size_t uOffset = vuCoefOffsets[uIndex];
						for(size_t bi = 0; bi < vusCoefCounts[uIndex]; bi++)
						{
							vdRefH[vusCoefBins[uOffset+bi]] += (WaveletSAT::typeSum)vdCoefSums[uOffset+bi];
						}
					}
				}

				// ADD-BY-LEETEN 2013/07/03-BEGIN
				if( iIsTestingBruteForce )
					continue;
				// ADD-BY-LEETEN 2013/07/03-END

				// truncate the numerical error
				double dError = 0.0;
				for(size_t b = 0; b < uNrOfBins; b++)
				{
					double dD = vdRefH[b] - vdH[b];
					dError += dD * dD;
				}
				double dRMSE = sqrt(dError / (double)uNrOfBins);

				if( dRMSE > 0.0 || iIsVerbose )
				{
					printf("Pos:");
					for(size_t d = 0; d < uNrOfDims; d++)
						  printf("%d,", vuBase[d]);
					printf("\n");

					for(size_t b = 0; b < uNrOfBins; b++)
					{
						if( vdH[b] == vdRefH[b] )
							printf( "\t\tH[%d]:%f\n", (unsigned int)b, vdH[b]);
						else
							printf( "\t\tH[%d]:%f,\t%f\n", (unsigned int)b, vdH[b], vdRefH[b]);
					}
					printf("RMSE:%f\n", dRMSE);
				}
			}

			LIBCLOCK_END(bIsPrintingTiming);
			long lAccumNrOfIORequests;		cSimpleNDFile._GetInteger(cSimpleNDFile.ACCUM_NR_OF_IO_REQUESTS,	&lAccumNrOfIORequests);	LOG_VAR(lAccumNrOfIORequests);
			long lMaxNrOfIORequests;		cSimpleNDFile._GetInteger(cSimpleNDFile.MAX_NR_OF_IO_REQUESTS,		&lMaxNrOfIORequests);	LOG_VAR(lMaxNrOfIORequests);
			long lMinNrOfIORequests;		cSimpleNDFile._GetInteger(cSimpleNDFile.MIN_NR_OF_IO_REQUESTS,		&lMinNrOfIORequests);	LOG_VAR(lMinNrOfIORequests);
		}
		else
		// ADD-BY-LEETEN 03/28/2013-END
		if(!szVolFilePath)
		{
			LIBCLOCK_BEGIN(bIsPrintingTiming);	// ADD-BY-LEETEN 12/30/2012

			cSimpleNDFile._SetInteger(cSimpleNDFile.RESET_IO_COUNTERS, 0);
			////////////////////////////////////////////////////////////////////////////
			// Now we can start to query SATs
			// load the data for testing

			size_t uNrOfDims = cSimpleNDFile.UGetNrOfDims();
			size_t uNrOfTestingValues = (size_t)iNrOfTestingValues;
			size_t uWinSize = iQueryWinLength;
			size_t uNrOfBins = cSimpleNDFile.UGetNrOfBins();	// it will be decided later

			vector<size_t> vuDimLengths;
			cSimpleNDFile._GetDataSize(vuDimLengths);
			size_t uNrOfValues = 1;	// ADD-BY-LEETEN 09/07/2012
			for(size_t d = 0; d < uNrOfDims; d++)
				uNrOfValues *= vuDimLengths[d];
			LIBCLOCK_END(bIsPrintingTiming);

			LIBCLOCK_BEGIN(bIsPrintingTiming);
			for(size_t t = 0; t < uNrOfTestingValues; t++)
			{
				vector<size_t> vuBase;
				vuBase.resize(uNrOfDims);
				for(size_t d = 0; d < uNrOfDims; d++)
					vuBase[d] = uWinSize + rand() % (vuDimLengths[d] - uWinSize);

				vector<WaveletSAT::typeSum> vdH;
				vdH.resize(uNrOfBins);

				vector<size_t> vuOffset;
				vuOffset.resize(uNrOfDims);
				for(size_t d = 0; d < uNrOfDims; d++)
					vuOffset[d] = vuBase[d] - uWinSize;
				cSimpleNDFile._GetRegionSums(vuOffset, vuBase, vdH);
			}

			LIBCLOCK_END(bIsPrintingTiming);
			long lAccumNrOfIORequests;		cSimpleNDFile._GetInteger(cSimpleNDFile.ACCUM_NR_OF_IO_REQUESTS,	&lAccumNrOfIORequests);	LOG_VAR(lAccumNrOfIORequests);
			long lMaxNrOfIORequests;		cSimpleNDFile._GetInteger(cSimpleNDFile.MAX_NR_OF_IO_REQUESTS,		&lMaxNrOfIORequests);	LOG_VAR(lMaxNrOfIORequests);
			long lMinNrOfIORequests;		cSimpleNDFile._GetInteger(cSimpleNDFile.MIN_NR_OF_IO_REQUESTS,		&lMinNrOfIORequests);	LOG_VAR(lMinNrOfIORequests);
		}
		else
		{
		LIBCLOCK_BEGIN(bIsPrintingTiming);	// ADD-BY-LEETEN 12/30/2012

		cSimpleNDFile._SetInteger(cSimpleNDFile.RESET_IO_COUNTERS, 0);
		////////////////////////////////////////////////////////////////////////////
		// Now we can start to query SATs
		// load the data for testing
		LOG_VAR(szVolFilePath);
		_ReadVolume(szVolFilePath);

		size_t uNrOfDims = (size_t)nin->dim;
		size_t uNrOfTestingValues = (size_t)iNrOfTestingValues;
		size_t uWinSize = iQueryWinLength;
		size_t uNrOfBins = cSimpleNDFile.UGetNrOfBins();	// it will be decided later

		vector<size_t> vuDimLengths;
		size_t uNrOfValues = 1;	// ADD-BY-LEETEN 09/07/2012
		for(size_t d = 0; d < uNrOfDims; d++)
		{
			size_t uDimLength = (size_t)nin->axis[d].size;
			vuDimLengths.push_back( uDimLength );
			uNrOfValues *= uDimLength;
		}

		cSimpleNDFile._SetHistogram(dValueMin, dValueMax);

		// decide the threshld to filter numerical error
		double dThreshold = cSimpleNDFile.DGetThreshold();

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
			vector< pair<WaveletSAT::typeBin, WaveletSAT::typeSum> > vuBins;
			cSimpleNDFile._MapValueToBins(vuBase, vdData[uIndex], vuBins);

			size_t uValueBin = vuBins[0].first;
			if( iIsVerbose )
				printf(")=\t%d,\n", (int)uValueBin);	// vdData[uIndex]);	// 

			vector<WaveletSAT::typeSum> vdH;

			vdH.resize(uNrOfBins);

			vector<size_t> vuOffset;
			vuOffset.resize(uNrOfDims);
			for(size_t d = 0; d < uNrOfDims; d++)
				vuOffset[d] = vuBase[d] - uWinSize;
			cSimpleNDFile._GetRegionSums(vuOffset, vuBase, vdH);

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
		long lAccumNrOfIORequests;		cSimpleNDFile._GetInteger(cSimpleNDFile.ACCUM_NR_OF_IO_REQUESTS,	&lAccumNrOfIORequests);	LOG_VAR(lAccumNrOfIORequests);
		long lMaxNrOfIORequests;		cSimpleNDFile._GetInteger(cSimpleNDFile.MAX_NR_OF_IO_REQUESTS,		&lMaxNrOfIORequests);	LOG_VAR(lMaxNrOfIORequests);
		long lMinNrOfIORequests;		cSimpleNDFile._GetInteger(cSimpleNDFile.MIN_NR_OF_IO_REQUESTS,		&lMinNrOfIORequests);	LOG_VAR(lMinNrOfIORequests);
	}
	}	// ADD-By-LEETEN 02/19/2013

	// ADD-BY-LEETEN 2013/12/01-BEGIN
    if( iIsComputingEntropy )
    {
		ASSERT_OR_LOG(szEntropyFilepathPrefix, "");

		/////////////////////////////////////////////////////////////
		LIBCLOCK_BEGIN(bIsPrintingTiming);
		cSimpleNDFile._SetInteger(cSimpleNDFile.PRINT_DECODE_BIN_TIMING, bIsPrintingTiming);
		vector<int> viLeft, viRight;
		for(size_t d = 0; d < cSimpleNDFile.UGetNrOfDims(); d++)
		{
			viLeft.push_back(-iEntropyWinRadius);
			viRight.push_back(+iEntropyWinRadius);
		}
		vector<double> valEntropyField;
#if	1	// TMP-MOD
		cSimpleNDFile._ComputeEntropy(viLeft, viRight, iNrOfEntropyBins, valEntropyField);
#else
		cSimpleNDFile._ComputeMedian(viLeft, viRight, iNrOfEntropyBins, valEntropyField);
#endif
		LIBCLOCK_END(bIsPrintingTiming);	

		/////////////////////////////////////////////////////////////
		vector<size_t> vuDimLengths;
		cSimpleNDFile._GetDataSize(vuDimLengths);

		size_t uDataSize = 1;
		for(size_t d = 0; d < vuDimLengths.size(); d++)
			uDataSize *= vuDimLengths[d];
		vector<float> vfNormalizedEntropy;
		vfNormalizedEntropy.resize(uDataSize);

		for(size_t i = 0; i < uDataSize; i++) 
			vfNormalizedEntropy[i] = (float)valEntropyField[i];

		WaveletSAT::_SaveNrrd<float>(
			vuDimLengths,
			vfNormalizedEntropy.data(),
			szEntropyFilepathPrefix
			);
	}
	// ADD-BY-LEETEN 2013/12/01-END

	LIBCLOCK_BEGIN(bIsPrintingTiming);
	cSimpleNDFile._ShowStatistics();
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_PRINT(bIsPrintingTiming);
	return 0;
}
