#include<vector>
using namespace std;
#include <stdio.h>
#include <assert.h>
#include "libclock.h"
#include "liblog.h"

template<typename T>
T
_DotProduct
(
	const vector<T>& v1,
	const vector<T>& v2
)
{
	T _Result(0);
	size_t uNrOfElements = min(v1.size(), v2.size());
	for(size_t e = 0; e < uNrOfElements; e++)
	{
		_Result += v1[e] * v2[e];
	}
	return _Result;
}

template<typename T>
T
_Sum
(
	const vector<T>& v1,
	const vector<T>& v2
)
{
	T _Result(0);
	size_t uNrOfElements = min(v1.size(), v2.size());
	for(size_t e = 0; e < uNrOfElements; e++)
	{
		_Result += v1[e] + v2[e];
	}
	return _Result;
}

template<typename T1, typename T2>
void
_Copy(vector<T1>& v1, const vector<T2>& v2)
{
	v1.clear();
	for(size_t e = 0; e < v2.size(); e++)
		v1.push_back((T1)v2[e]);
}


int
main(int argn, char* argv[])
{
	size_t uNrOfElements = 100000;
	size_t uNrOfTests = 10000;

	LIBCLOCK_INIT(1, "main");

	// test int
	vector<int> vi1;
	vi1.resize(uNrOfElements);
	vector<int> vi2;
	vi2.resize(uNrOfElements);
	for(size_t e = 0; e < uNrOfElements; e++)
	{
		vi1[e] = rand() % 128;
		vi2[e] = rand() % 128;
	}
	LIBCLOCK_BEGIN(1);	
	int iProd;
	for(size_t t = 0; t < uNrOfTests; t++)
		iProd = _DotProduct<int>(vi1, vi2);	
	LIBCLOCK_END(1);
	LOG_VAR(iProd);

	int iSum;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		iSum = _Sum<int>(vi1, vi2);	
	LIBCLOCK_END(1);
	LOG_VAR(iSum);

	// test unsigned int
	vector<unsigned int> vui1;		_Copy<unsigned int, int>(vui1, vi1);
	vector<unsigned int> vui2;		_Copy<unsigned int, int>(vui2, vi2);
	unsigned int uiProd;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		uiProd = _DotProduct<unsigned int>(vui1, vui2);	
	LIBCLOCK_END(1);
	LOG_VAR(uiProd);

	unsigned int uiSum;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		uiSum = _Sum<unsigned int>(vui1, vui2);	
	LIBCLOCK_END(1);
	LOG_VAR(uiSum);

	// test long
	vector<long> vl1;		_Copy<long, int>(vl1, vi1);
	vector<long> vl2;		_Copy<long, int>(vl2, vi2);
	long lProd;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		lProd = _DotProduct<long>(vl1, vl2);	
	LIBCLOCK_END(1);
	LOG_VAR(lProd);

	long lSum;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		lSum = _Sum<long>(vl1, vl2);	
	LIBCLOCK_END(1);
	LOG_VAR(lSum);

	// test long
	vector<unsigned long> vul1;		_Copy<unsigned long, int>(vul1, vi1);
	vector<unsigned long> vul2;		_Copy<unsigned long, int>(vul2, vi2);
	long ulProd;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		ulProd = _DotProduct<unsigned long>(vul1, vul2);	
	LIBCLOCK_END(1);
	LOG_VAR(ulProd);

	long ulSum;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		ulSum = _Sum<unsigned long>(vul1, vul2);	
	LIBCLOCK_END(1);
	LOG_VAR(ulSum);

	// test size_t
	vector<size_t> vu1;		_Copy<size_t, int>(vu1, vi1);
	vector<size_t> vu2;		_Copy<size_t, int>(vu2, vi2);
	size_t uProd;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		uProd = _DotProduct<size_t>(vu1, vu2);	
	LIBCLOCK_END(1);
	LOG_VAR(uProd);

	size_t uSum;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		uSum = _Sum<size_t>(vu1, vu2);	
	LIBCLOCK_END(1);
	LOG_VAR(uSum);

	// test float
	vector<float> vf1;		_Copy<float, int>(vf1, vi1);
	vector<float> vf2;		_Copy<float, int>(vf2, vi2);
	float fProd;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		fProd = _DotProduct<float>(vf1, vf2);	
	LIBCLOCK_END(1);
	LOG_VAR(fProd);

	float fSum;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		fSum = _Sum<float>(vf1, vf2);	
	LIBCLOCK_END(1);
	LOG_VAR(fSum);

	// test double
	vector<double> vd1;		_Copy<double, int>(vd1, vi1);
	vector<double> vd2;		_Copy<double, int>(vd2, vi2);
	double dProd;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		dProd = _DotProduct<double>(vd1, vd2);	
	LIBCLOCK_END(1);
	LOG_VAR(dProd);

	double dSum;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		dSum = _Sum<double>(vd1, vd2);	
	LIBCLOCK_END(1);
	LOG_VAR(dSum);

	LIBCLOCK_PRINT(1);

	return 0;
}