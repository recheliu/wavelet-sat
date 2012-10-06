#include <math.h>
#include <assert.h>
#include <stdlib.h> // ADD-BY-LEETEN 09/09/2012
#include "libclock.h"
#include "libopt.h"

#include "SimpleND.h"

int
main(int argn, char* argv[])
{
	bool bIsPrintingTiming = true;
	LIBCLOCK_INIT(bIsPrintingTiming, __FUNCTION__);
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	CSimpleND<int> cSimpleND;

	_OPTInit();			// initialize the option parser

	enum {
		DATA_SAMPLING_BLOCK,
		DATA_SAMPLING_RANDOM_BLOCK,	// ADD-BY-LEETEN 09/14/2012
		DATA_SAMPLING_RANDOM,
		DATA_SAMPLING_JUMP,
		NR_OF_DATA_SAMPLING,
		DATA_SAMPLING_DEFAULT = DATA_SAMPLING_BLOCK
	};
	int iDataSampling;
	_OPTAddEnum("--data-sampling", &iDataSampling, DATA_SAMPLING_DEFAULT, NR_OF_DATA_SAMPLING,
		"block",	DATA_SAMPLING_BLOCK,
		"random-block",	DATA_SAMPLING_RANDOM_BLOCK,	// ADD-BY-LEETEN 09/14/2012
		"random",	DATA_SAMPLING_RANDOM,
		"jump",		DATA_SAMPLING_JUMP,
		NULL);

	int iDimLength = 128;
	_OPTAddIntegerVector(
		"--dim-length", 1,
		&iDimLength, iDimLength);

	int iNrOfDims = 2;
	_OPTAddIntegerVector(
		"--n-dims", 1,
		&iNrOfDims, iNrOfDims);

	int iMaxLevel = 0;
	_OPTAddIntegerVector(
		"--max-level", 1,
		&iMaxLevel, iMaxLevel);

	int iValueMax = 32;
	_OPTAddIntegerVector(
		"--value-max", 1,
		&iValueMax, iValueMax);

	int iNrOfBins = 8;	// iValueMax;
	_OPTAddIntegerVector(
		"--n-bins", 1,
		&iNrOfBins, iNrOfBins);

	int iIsTestingQuery = 0; 
	_OPTAddBoolean(
		"--is-testing-query", &iIsTestingQuery, iIsTestingQuery);

	int iIsVerbose = 0; 
	_OPTAddBoolean(
		"--is-verbose", &iIsVerbose, iIsVerbose);

	// ADD-BY-LEETEN 09/14/2012-BEGIN
	int iSizeOfFullArrays = 0;
	_OPTAddIntegerVector(
		"--size-of-full-arrays", 1,
		&iSizeOfFullArrays, iSizeOfFullArrays);
	_OPTAddComment("--size-of-full-arrays", 
		"Size (in MB) of the full arrays from all bin SATs");
	// ADD-BY-LEETEN 09/14/2012-END

	bool bIsOptParsed = BOPTParse(argv, argn, 1);
	assert(bIsOptParsed);

	size_t uDimLength = iDimLength;
	size_t uNrOfDims = iNrOfDims;
	size_t uMaxLevel = iMaxLevel;
	size_t uWinSize = 1;
	size_t uNrOfBins = iNrOfBins;
	// ADD-By-LEETEN 09/07/2012-END

	// ADD-BY-LEETEN 09/14/2012-BEGIN
	cSimpleND._SetLong(CSimpleND<int>::SIZE_OF_FULL_ARRAYS, (long)iSizeOfFullArrays);

	// create a lookup table to shuffule the value
	vector<int> viShuffleTable;
	for(size_t i = 0; i < iValueMax; i++)
		viShuffleTable.push_back((int)i);

	for(size_t i = iValueMax; i > 0; i--)
	{
		size_t uRand = rand() % i;
		swap(viShuffleTable[i - 1], viShuffleTable[uRand]);
	}
	// ADD-BY-LEETEN 09/14/2012-END

	size_t uNrOfTestingValues = uDimLength;

	// Step 1: Setup up the data size
	vector<size_t> vuDimLengths;
	size_t uNrOfValues = 1;	// ADD-BY-LEETEN 09/07/2012
	for(size_t d = 0; d < uNrOfDims; d++)
	{
		vuDimLengths.push_back(uDimLength);
		uNrOfValues *= uDimLength;
	}
	cSimpleND._SetDimLengths(vuDimLengths);

	// Step 2: Allocate the needed #SATs
	cSimpleND._SetHistogram(uNrOfBins, 0, iValueMax);
	vector<size_t> vuDimMaxLevels;
	for(size_t d = 0; d < uNrOfDims; d++)
		vuDimMaxLevels.push_back(uMaxLevel);
	cSimpleND._AllocateBins(uNrOfBins, vuDimMaxLevels);
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_BEGIN(bIsPrintingTiming);
	// Step 3: Add the value to the SAT
	vector<size_t> vuValueBins;
	for(size_t i = 0; i < uNrOfValues; i++)
	{
		vector<size_t> vuPos;
		for(size_t 
			d = 0, uCoord = i; 
			d < uNrOfDims; 
			uCoord /= vuDimLengths[d], d++)
			vuPos.push_back(uCoord % vuDimLengths[d]);

		int iValue = 0;
		switch(iDataSampling)
		{
		case DATA_SAMPLING_RANDOM:
			iValue = rand()%iValueMax;
			break;

		case DATA_SAMPLING_RANDOM_BLOCK:	// ADD-BY-LEETEN 09/14/2012
		case DATA_SAMPLING_BLOCK:
			{
			int iValueDimMax = (int)floor(pow((double)iValueMax, 1.0/(double)uNrOfDims));
			for(size_t 
				d = 0, uPrevDimSize = 1; 
				d < uNrOfDims; 
				uPrevDimSize *= iValueDimMax, d++)
				iValue += (int)uPrevDimSize * (int)floorf((float)iValueDimMax * (float)vuPos[d] / (float)vuDimLengths[d]);

			// ADD-BY-LEETEN 09/14/2012-BEGIN
			switch(iDataSampling)
			{
			case DATA_SAMPLING_RANDOM_BLOCK:
				iValue = (int)viShuffleTable[iValue];
				break;
			}
			// ADD-BY-LEETEN 09/14/2012-END
			}
			break;

		case DATA_SAMPLING_JUMP:
			{
			int iDimBin = (int)floor(pow((double)uNrOfBins, 1.0/(double)uNrOfDims));
			for(size_t 
				d = 0, uPrevValue = 1; 
				d < uNrOfDims; 
				uPrevValue *= iDimBin, d++)
				iValue += (int)uPrevValue * (int)(vuPos[d] % iDimBin);
			}
			break;
		}

		cSimpleND._AddValue(vuPos, iValue);

		vector< pair<size_t, double> > vuBins;
		cSimpleND._MapValueToBins(vuPos, iValue, vuBins);
		vuValueBins.push_back(vuBins[0].first);
		//		printf("%d\n", iValue);
	}

	// ADD-BY-LEETEN 09/07/2012-BEGIN
	// Step 4: Finalize the SAT computation
	cSimpleND._Finalize();
	// ADD-BY-LEETEN 09/07/2012-END

	LIBCLOCK_END(bIsPrintingTiming);

	////////////////////////////////////////////////////////////////////////////
	// Now we can start to query SATs
	LIBCLOCK_BEGIN(bIsPrintingTiming);

	cSimpleND._ShowStatistics();
	LIBCLOCK_END(bIsPrintingTiming);

	// ADD-BY-LEETEN	09/09/2012-BEGIN
	if(iIsTestingQuery)
	{
	// ADD-BY-LEETEN	09/09/2012-END
	double dThreshold = cSimpleND.DGetThreshold();	// ADD-BY-LEETEN 10/05/2012
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
	// ADD-BY-LEETEN 09/07/2012-BEGIN
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_BEGIN(bIsPrintingTiming);
	// ADD-BY-LEETEN 09/07/2012-END

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
		if( iIsVerbose )
		printf(")=%d,\n", (int)vuValueBins[uIndex]);

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

			// ADD-BY-LEETEN 09/09/2012-BEGIN
			#if	0	// DEL-BY-LEETEN 10/05/2012-BEGIN
				if( iIsVerbose )
					printf("\t%+d,", iSign);
			#endif		// DEL-BY-LEETEN 10/05/2012-END
			// ADD-BY-LEETEN 09/09/2012-END
			for(size_t b = 0; b < uNrOfBins; b++)
			{	// ADD-BY-LEETEN 09/09/2012
				vdH[b] += iSign * vdIH[b]; 
			// ADD-BY-LEETEN 09/09/2012-BEGIN

				#if	0	// DEL-BY-LEETEN 10/05/2012-BEGIN
					if( iIsVerbose )	// ADD-BY-LEETEN 09/12/2012
					printf( "%+.2f,", vdIH[b]);
				#endif		// DEL-BY-LEETEN 10/05/2012-END
			}
			#if	0	// DEL-BY-LEETEN 10/05/2012-BEGIN
				if( iIsVerbose )	
					printf("\n");
			#endif		// DEL-BY-LEETEN 10/05/2012-END
			// ADD-BY-LEETEN 09/09/2012-END
		}

		double dError = 0.0;
		#if	0	// DEL-BY-LEETEN 10/05/2012-BEGIN
			if( iIsVerbose )
			printf("H:");
		#endif		// DEL-BY-LEETEN 10/05/2012-END
		for(size_t b = 0; b < uNrOfBins; b++)
		{
			#if	0	// DEL-BY-LEETEN 10/05/2012-BEGIN
				if( iIsVerbose )
				printf( "%+.2f,", vdH[b]);
			#endif		// DEL-BY-LEETEN 10/05/2012-END

			if(b == vuValueBins[uIndex])
				dError += pow(1.0 - vdH[b], 2.0);
			else
				dError += pow(vdH[b], 2.0);
		}
		#if	0	// DEL-BY-LEETEN 10/05/2012-BEGIN
			if( iIsVerbose )	
			printf("E:%f\n", dError);
		#endif		// DEL-BY-LEETEN 10/05/2012-END

		// ADD-BY-LEETEN 10/05/2012-BEGIN
		if( iIsVerbose )
		{
			printf("H:");
			for(size_t b = 0; b < uNrOfBins; b++)
				if( fabs(vdH[b]) > dThreshold )
				  printf( "\t%d:%+.2f\n", (unsigned int)b, vdH[b]); // MOD-BY-LEETEN 10/02/2012-FROM: printf( "\t%d:%+.2f\n", b, vdH[b]);
			printf("E:%f\n", dError);
		}
		// ADD-BY-LEETEN 10/05/2012-END
	}
	LIBCLOCK_END(bIsPrintingTiming);
	}	// ADD-BY-LEETEN	09/09/2012-END

	LIBCLOCK_PRINT(bIsPrintingTiming);
	return 0;
}
