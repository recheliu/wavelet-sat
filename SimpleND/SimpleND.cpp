#include <math.h>
#include <assert.h>
#include <stdlib.h>
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
		DATA_SAMPLING_RANDOM_BLOCK,	
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

	int iSizeOfFullArrays = 0;
	_OPTAddIntegerVector(
		"--size-of-full-arrays", 1,
		&iSizeOfFullArrays, iSizeOfFullArrays);
	_OPTAddComment("--size-of-full-arrays", 
		"Size (in MB) of the full arrays from all bin SATs");

	bool bIsOptParsed = BOPTParse(argv, argn, 1);
	assert(bIsOptParsed);

	size_t uDimLength = iDimLength;
	size_t uNrOfDims = iNrOfDims;
	size_t uMaxLevel = iMaxLevel;
	size_t uWinSize = 1;
	size_t uNrOfBins = iNrOfBins;

	cSimpleND._SetLong(CSimpleND<int>::SIZE_OF_FULL_ARRAYS, (long)iSizeOfFullArrays);

	// create a lookup table to shuffule the value
	vector<int> viShuffleTable;
	// MOD-BY-LEETEN 10/17/2012-FROM:	for(size_t i = 0; i < iValueMax; i++)
	for(int i = 0; i < iValueMax; i++)
	// MOD-BY-LEETEN 10/17/2012-END
		viShuffleTable.push_back((int)i);

	// MOD-BY-LEETEN 10/17/2012-FROM:	for(size_t i = iValueMax; i > 0; i--)
	for(int i = iValueMax; i > 0; i--)
	// MOD-BY-LEETEN 10/17/2012-END
	{
		size_t uRand = rand() % i;
		swap(viShuffleTable[i - 1], viShuffleTable[uRand]);
	}

	size_t uNrOfTestingValues = uDimLength;

	// Step 1: Setup up the data size
	vector<size_t> vuDimLengths;
	size_t uNrOfValues = 1;	
	for(size_t d = 0; d < uNrOfDims; d++)
	{
		vuDimLengths.push_back(uDimLength);
		uNrOfValues *= uDimLength;
	}
	// MOD-BY-LEETEN 10/27/2012-FROM:	cSimpleND._SetDimLengths(vuDimLengths);
	cSimpleND._Set(vuDimLengths, uNrOfBins);
	// MOD-BY-LEETEN 10/27/2012-END

	// Step 2: Allocate the needed #SATs
	cSimpleND._SetHistogram(uNrOfBins, 0, iValueMax);
	#if	0	// MOD-BY-LEETEN 10/27/2012-FROM:
	vector<size_t> vuDimMaxLevels;
	for(size_t d = 0; d < uNrOfDims; d++)
		vuDimMaxLevels.push_back(uMaxLevel);
	cSimpleND._AllocateBins(uNrOfBins, vuDimMaxLevels);
	#else		// MOD-BY-LEETEN 10/27/2012-TO:
	cSimpleND._Allocate();
	#endif		// MOD-BY-LEETEN 10/27/2012-END
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

		case DATA_SAMPLING_RANDOM_BLOCK:
		case DATA_SAMPLING_BLOCK:
			{
			int iValueDimMax = (int)floor(pow((double)iValueMax, 1.0/(double)uNrOfDims));
			for(size_t 
				d = 0, uPrevDimSize = 1; 
				d < uNrOfDims; 
				uPrevDimSize *= iValueDimMax, d++)
				iValue += (int)uPrevDimSize * (int)floorf((float)iValueDimMax * (float)vuPos[d] / (float)vuDimLengths[d]);

			switch(iDataSampling)
			{
			case DATA_SAMPLING_RANDOM_BLOCK:
				iValue = (int)viShuffleTable[iValue];
				break;
			}
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

	// Step 4: Finalize the SAT computation
	cSimpleND._Finalize();

	LIBCLOCK_END(bIsPrintingTiming);

	////////////////////////////////////////////////////////////////////////////
	// Now we can start to query SATs
	LIBCLOCK_BEGIN(bIsPrintingTiming);

	cSimpleND._ShowStatistics();
	LIBCLOCK_END(bIsPrintingTiming);

	if(iIsTestingQuery)
	{
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
		if( iIsVerbose )
			printf(")=\t%d,\n", (int)vuValueBins[uIndex]);	

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
			{	
				vdH[b] += iSign * vdIH[b]; 
			}
		}

		double dError = 0.0;
		for(size_t b = 0; b < uNrOfBins; b++)
		{
			if(b == vuValueBins[uIndex])
				dError += pow(1.0 - vdH[b], 2.0);
			else
				dError += pow(vdH[b], 2.0);
		}

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
