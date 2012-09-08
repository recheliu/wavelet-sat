#include <math.h>
#include <assert.h>
#include "libclock.h"

#include "WaveletSAT.h"

#define PRINT_OUTPUT	0	// ADD-BY-LEETEN 09/07/2012

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
		vector<pair<size_t, double>>& vpBins,
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

int
main(int arhn, char* argv[])
{
	bool bIsPrintingTiming = true;
	LIBCLOCK_INIT(bIsPrintingTiming, __FUNCTION__);
	LIBCLOCK_BEGIN(bIsPrintingTiming);
	CSimpleND<int> cSimpleND;

	size_t uDimLength = 128;
	size_t uNrOfDims = 2;
	// DEL-BY-LEETEN 09/07/2012:	size_t uNrOfValues = 1;
	size_t uMaxLevel = 0;
	size_t uWinSize = 1;
	// ADD-By-LEETEN 09/07/2012-BEGIN
	int iValueMax = 32;
	size_t uNrOfBins = 8;	// iValueMax;
	// ADD-By-LEETEN 09/07/2012-END

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

		int iValue = rand()%iValueMax;
		cSimpleND._AddValue(vuPos, iValue);

		#if	0	// MOD-BY-LEETEN 09/07/2012-FROM:
			vector<size_t> vuBins;
			cSimpleND._MapValueToBins(vuPos, iValue, vuBins);
			vuValueBins.push_back(vuBins[0]);
		#else		// MOD-BY-LEETEN 09/07/2012-TO:
		vector<pair<size_t, double>> vuBins;
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

	vector<double> vdBinEnergies;
	cSimpleND._GetEnergy
	(
		vdBinEnergies
	);

	#if	PRINT_OUTPUT		// ADD-BY-LEETEN 09/07/2012
	for(size_t b = 0; b < vdBinEnergies.size(); b++)
		printf("Bin (%d): %f\n", b, vdBinEnergies[b]);
	#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_BEGIN(bIsPrintingTiming);

	size_t uNrOfIHs = 1 << uNrOfDims;
	vector<vector<size_t>> vvuOffsets;
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
		#if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
		printf("B(");
		#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
		for(size_t 
			d = 0, uDimLengthProduct = 1; 
			d < uNrOfDims; 
			uDimLengthProduct *= vuDimLengths[d], d++)
		{
			size_t uPos = uWinSize + rand() % (vuDimLengths[d] - uWinSize);
			vuBase.push_back(uPos);
			uIndex += uPos * uDimLengthProduct;

			#if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
			printf("%4d,", uPos);
			#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
		}
		#if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
		printf(")=%d,", vuValueBins[uIndex]);
		#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012

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
		#if	PRINT_OUTPUT		// ADD-BY-LEETEN 09/07/2012
		printf("H:");
		#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
		for(size_t b = 0; b < uNrOfBins; b++)
		{
			#if	PRINT_OUTPUT		// ADD-BY-LEETEN 09/07/2012
			printf( "%+.2f,", vdH[b]);
			#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
			if(b == vuValueBins[uIndex])
				dError += pow(1.0 - vdH[b], 2.0);
			else
				dError += pow(vdH[b], 2.0);
		}
		#if	PRINT_OUTPUT		// ADD-BY-LEETEN 09/07/2012
		printf("E:%f\n", dError);
		#endif	// #if	PRINT_OUTPUT	// ADD-BY-LEETEN 09/07/2012
	}
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_PRINT(bIsPrintingTiming);
	return 0;
}