#include <math.h>
#include <assert.h>
#include "libclock.h"

#include "WaveletSAT.h"

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
		vector<size_t>& vuBins,
		void *_Reserved = NULL
	)
	{
		T clampedValue = min(max(value, valueMin), valueMax);
		size_t uBin = (size_t)floorf((float)(uNrOfBins * (clampedValue - valueMin))/(float)(valueMax - valueMin));
		uBin = min(uBin, uNrOfBins - 1);
		vuBins.clear();
		vuBins.push_back(uBin);
	}
	
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
	size_t uNrOfValues = 1;
	size_t uMaxLevel = 0;
	size_t uWinSize = 1;

	size_t uNrOfTestingValues = uDimLength;

	int iValueMax = 32;
	size_t uNrOfBins = 8;	// iValueMax;

	vector<size_t> vuDimLengths;
	uNrOfValues = 1;
	for(size_t d = 0; d < uNrOfDims; d++)
	{
		vuDimLengths.push_back(uDimLength);
		uNrOfValues *= uDimLength;
	}
	cSimpleND._SetDimLengths(vuDimLengths);

	cSimpleND._SetHistogram(uNrOfBins, 0, iValueMax);
	vector<size_t> vuDimMaxLevels;
	for(size_t d = 0; d < uNrOfDims; d++)
		vuDimMaxLevels.push_back(uMaxLevel);
	cSimpleND._AllocateBins(uNrOfBins, vuDimMaxLevels);
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_BEGIN(bIsPrintingTiming);
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

		vector<size_t> vuBins;
		cSimpleND._MapValueToBins(vuPos, iValue, vuBins);
		vuValueBins.push_back(vuBins[0]);
		//		printf("%d\n", iValue);
	}
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_BEGIN(bIsPrintingTiming);

	vector<double> vdBinEnergies;
	cSimpleND._GetEnergy
	(
		vdBinEnergies
	);

	for(size_t b = 0; b < vdBinEnergies.size(); b++)
		printf("Bin (%d): %f\n", b, vdBinEnergies[b]);
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

	for(size_t t = 0; t < uNrOfTestingValues; t++)
	{
		vector<size_t> vuBase;
		size_t uIndex = 0;
		printf("B(");
		for(size_t 
			d = 0, uDimLengthProduct = 1; 
			d < uNrOfDims; 
			uDimLengthProduct *= vuDimLengths[d], d++)
		{
			size_t uPos = uWinSize + rand() % (vuDimLengths[d] - uWinSize);
			vuBase.push_back(uPos);
			uIndex += uPos * uDimLengthProduct;

			printf("%4d,", uPos);
		}
		printf(")=%d,", vuValueBins[uIndex]);

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
		printf("H:");
		for(size_t b = 0; b < uNrOfBins; b++)
		{
			printf( "%+.2f,", vdH[b]);
			if(b == vuValueBins[uIndex])
				dError += pow(1.0 - vdH[b], 2.0);
			else
				dError += pow(vdH[b], 2.0);
		}
		printf("E:%f\n", dError);
	}
	LIBCLOCK_END(bIsPrintingTiming);

	LIBCLOCK_PRINT(bIsPrintingTiming);
	return 0;
}