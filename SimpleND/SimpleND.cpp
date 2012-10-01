#include <math.h>
#include <assert.h>
#include <stdlib.h> // ADD-BY-LEETEN 09/09/2012
#include "libclock.h"
#include "libopt.h"

#if	0	// MOD-BY-LEETEN 09/30/2012-FROM:
#include "WaveletSAT.h"

// DEL-BY-LEETEN	09/09/2012:	#define PRINT_OUTPUT	1	// ADD-BY-LEETEN 09/07/2012
using namespace WaveletSAT;

template<class T>
class CSimpleND:
	public CBase<T>
{
	size_t uNrOfBins;
	T valueMin, valueMax;
public:
	virtual 
	void 
	_MapValueToBins
	(
		const vector<size_t>& vuPos,
		const T& value, 
		// MOD-BY-LEETEN 09/07/2012-FROM:	vector<size_t>& vuBins,
		// MOD-BY-LEETEN 09/09/2012-FROM:	vector<pair<size_t, double>>& vpBins,
		vector< pair<size_t, double> >& vpBins,
		// MOD-BY-LEETEN 09/09/2012-END
		// MOD-BY-LEETEN 09/07/2012-END
		void *_Reserved = NULL
	)
	{
		T clampedValue = min(max(value, valueMin), valueMax);
		size_t uBin = (size_t)floorf((float)(uNrOfBins * (clampedValue - valueMin))/(float)(valueMax - valueMin));
		uBin = min(uBin, uNrOfBins - 1);
		#if	0	// MOD-BY-LEETEN 09/07/2012-FORM:
			vuBins.clear();
			vuBins.push_back(uBin);
		#else		// MOD-BY-LEETEN 09/07/2012-TO:
		vpBins.clear();
		vpBins.push_back(pair<size_t, double>(uBin, 1.0));
		#endif		// MOD-BY-LEETEN 09/07/2012-END
	}
	
	#if	0	// DEL-BY-LEETEN 09/07/2012-BEGIN
		virtual
		void
		_MapValueToBinWeight
		(
			const vector<size_t>& vuPos,
			const T& value, 
			size_t uBin,
			double& dW,
			void *_Reserved = NULL
		)
		{
			dW = 1.0;
		}
	#endif	// DEL-BY-LEETEN 09/07/2012-END

	////////////////////////////////////////////////////////////
	void
	_SetHistogram
	(
		size_t uNrOfBins,
		const T& valueMin, 
		const T& valueMax
	)
	{
		this->uNrOfBins = uNrOfBins;
		this->valueMin = valueMin;
		this->valueMax = valueMax;
	}

	void
	_AddValue
	(
		const vector<size_t>& vuPos,
		const T& value,
		void *_Reserved = NULL
	)
	{
		_Update(vuPos, value);
	}
};
#else	// MOD-BY-LEETEN 09/30/2012-TO:

#include "SimpleND.h"

#endif	// MOD-BY-LEETEN 09/30/2012-END

int
// MOD-BY-LEETEN	09/09/2012-FROM:	main(int arhn, char* argv[])
main(int argn, char* argv[])
// DEL-BY-LEETEN	09/09/2012-END
{
	bool bIsPrintingTiming = true;
	LIBCLOCK_INIT(bIsPrintingTiming, __FUNCTION__);
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	CSimpleND<int> cSimpleND;

	#if	0	// MOD-BY-LEETEN	09/09/2012-FROM:
		size_t uDimLength = 128;
		size_t uNrOfDims = 2;
		// DEL-BY-LEETEN 09/07/2012:	size_t uNrOfValues = 1;
		size_t uMaxLevel = 0;
		size_t uWinSize = 1;
		// ADD-By-LEETEN 09/07/2012-BEGIN
		int iValueMax = 32;
		size_t uNrOfBins = 8;	// iValueMax;
		// ADD-By-LEETEN 09/07/2012-END
	#else		// MOD-BY-LEETEN	09/09/2012-TO:
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
	#endif	// MOD-BY-LEETEN	09/09/2012-END

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
	#if	0	// DEL-BY-LEETEN 09/07/2012-BEGIN
		int iValueMax = 32;
		size_t uNrOfBins = 8;	// iValueMax;
	#endif		// DEL-BY-LEETEN 09/07/2012-END

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

		// MOD-BY-LEETEN	09/09/2012-FROM:		int iValue = rand()%iValueMax;
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
		// MOD-BY-LEETEN	09/09/2012-END

		cSimpleND._AddValue(vuPos, iValue);

		#if	0	// MOD-BY-LEETEN 09/07/2012-FROM:
			vector<size_t> vuBins;
			cSimpleND._MapValueToBins(vuPos, iValue, vuBins);
			vuValueBins.push_back(vuBins[0]);
		#else		// MOD-BY-LEETEN 09/07/2012-TO:
		// MOD-BY-LEETEN 09/09/2012-FROM:		vector<pair<size_t, double>> vuBins;
		vector< pair<size_t, double> > vuBins;
		// MOD-BY-LEETEN 09/09/2012-END
		cSimpleND._MapValueToBins(vuPos, iValue, vuBins);
		vuValueBins.push_back(vuBins[0].first);
		#endif		// MOD-BY-LEETEN 09/07/2012-END
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

	#if	0	// MOD-BY-LEETEN	09/09/2012-FROM:
		vector<double> vdBinEnergies;
		cSimpleND._GetEnergy
		(
			vdBinEnergies
		);

		#if	PRINT_OUTPUT		// ADD-BY-LEETEN 09/07/2012
		for(size_t b = 0; b < vdBinEnergies.size(); b++)
			printf("Bin (%d): %f\n", b, vdBinEnergies[b]);
		#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
	#else		// MOD-BY-LEETEN	09/09/2012-TO:
	cSimpleND._ShowStatistics();
	#endif		// MOD-BY-LEETEN	09/09/2012-END
	LIBCLOCK_END(bIsPrintingTiming);

	// ADD-BY-LEETEN	09/09/2012-BEGIN
	if(iIsTestingQuery)
	{
	// ADD-BY-LEETEN	09/09/2012-END
	LIBCLOCK_BEGIN(bIsPrintingTiming);

	size_t uNrOfIHs = 1 << uNrOfDims;
	// MOD-BY-LEETEN 09/09/2012-FROM:	vector<vector<size_t>> vvuOffsets;
	vector< vector<size_t> > vvuOffsets;
	// MOD-BY-LEETEN 09/09/2012-END
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
		if( iIsVerbose )	// MOD-BY-LEETEN	09/09/2012-FROM:	#if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
		printf("B(");
		// DEL-BY-LEETEN	09/09/2012:	#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
		for(size_t 
			d = 0, uDimLengthProduct = 1; 
			d < uNrOfDims; 
			uDimLengthProduct *= vuDimLengths[d], d++)
		{
			size_t uPos = uWinSize + rand() % (vuDimLengths[d] - uWinSize);
			vuBase.push_back(uPos);
			uIndex += uPos * uDimLengthProduct;

			if( iIsVerbose )	// MOD-BY-LEETEN	09/09/2012-FROM:	#if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
#if 0 // MOD-BY-LEETEN 09/09/2012-FROM:
			  
			printf("%3d,", uPos);	// MOD-BY-LEETEN	09/09/2012-FROM:	printf("%4d,", uPos);
#else // MOD-BY-LEETEN 09/09/2012-TO:
			printf("%3d,", (int)uPos);
#endif // MOD-BY-LEETEN 09/09/2012-END
			// DEL-BY-LEETEN	09/09/2012:	#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
		}
		if( iIsVerbose )	// MOD-BY-LEETEN	09/09/2012-FROM:	#if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
		// MOD-BY-LEETEN 09/09/2012-FROM:		printf(")=%d,", vuValueBins[uIndex]);
		  // MOD-BY-LEETEN 09/09/2012-FROM:		printf(")=%d,\n", vuValueBins[uIndex]);
		printf(")=%d,\n", (int)vuValueBins[uIndex]);
		// MOD-BY-LEETEN 09/09/2012-END
		// MOD-BY-LEETEN 09/09/2012-END
		// DEL-BY-LEETEN	09/09/2012:	#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012

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
			if( iIsVerbose )
				printf("\t%+d,", iSign);
			// ADD-BY-LEETEN 09/09/2012-END
			for(size_t b = 0; b < uNrOfBins; b++)
			{	// ADD-BY-LEETEN 09/09/2012
				vdH[b] += iSign * vdIH[b]; 
			// ADD-BY-LEETEN 09/09/2012-BEGIN
				if( iIsVerbose )	// ADD-BY-LEETEN 09/12/2012
				printf( "%+.2f,", vdIH[b]);
			}
			if( iIsVerbose )	
				printf("\n");
			// ADD-BY-LEETEN 09/09/2012-END
		}

		double dError = 0.0;
		if( iIsVerbose )	// MOD-BY-LEETEN	09/09/2012-FROM:	#if	PRINT_OUTPUT		// ADD-BY-LEETEN 09/07/2012
		printf("H:");
		// DEL-BY-LEETEN	09/09/2012:	#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
		for(size_t b = 0; b < uNrOfBins; b++)
		{
			if( iIsVerbose )	// MOD-BY-LEETEN	09/09/2012-FROM:	#if	PRINT_OUTPUT		// ADD-BY-LEETEN 09/07/2012
			printf( "%+.2f,", vdH[b]);
			// DEL-BY-LEETEN	09/09/2012:	#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
			if(b == vuValueBins[uIndex])
				dError += pow(1.0 - vdH[b], 2.0);
			else
				dError += pow(vdH[b], 2.0);
		}
		if( iIsVerbose )	// MOD-BY-LEETEN	09/09/2012-FROM:	#if	PRINT_OUTPUT		// ADD-BY-LEETEN 09/07/2012
		printf("E:%f\n", dError);
		// DEL-BY-LEETEN	09/09/2012:	#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
	}
	LIBCLOCK_END(bIsPrintingTiming);
	}	// ADD-BY-LEETEN	09/09/2012-END

	LIBCLOCK_PRINT(bIsPrintingTiming);
	return 0;
}
