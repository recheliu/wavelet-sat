#include <typeinfo>	// ADD-BY-LEETEN 10/28/2012
#include<vector>
using namespace std;
#include <stdio.h>
#include <assert.h>
#include <stdlib.h> // ADD-BY-LEETEN 10/22/2012
#include "libclock.h"
#include "liblog.h"

// ADD-BY-LEETEN 10/28/2012-BEGIN
template<typename T>
T
_Multiply
(
	const vector<T>& v
)
{
	T _Result = v[0];
	for(size_t e = 1; e < v.size(); e++)
	{
		_Result *= v[e];
	}
	return _Result;
}
// ADD-BY-LEETEN 10/28/2012-END

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

// ADD-BY-LEETEN 10/28/2012-BEGIN
vector<int> vi1;
vector<int> vi2;

template<typename T>
void
_TestType()
{
	size_t uNrOfTests = 10000;

	LIBCLOCK_INIT(1, typeid(T).name());	
	LIBCLOCK_BEGIN(1);	
	vector<T> v1;
	v1.resize(vi1.size());
	_Copy<T, int>(v1, vi1);

	vector<T> v2;
	v2.resize(vi2.size());
	_Copy<T, int>(v2, vi2);
	LIBCLOCK_END(1);

	LIBCLOCK_BEGIN(1);	
	T Prod;
	for(size_t t = 0; t < uNrOfTests; t++)
		Prod = _DotProduct<T>(v1, v2);	
	LIBCLOCK_END(1);
	LOG_VAR(Prod);

	T Sum;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		Sum = _Sum<T>(v1, v2);	
	LIBCLOCK_END(1);
	LOG_VAR(Sum);

	T Mul;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		Mul = _Multiply<T>(v1);	
	LIBCLOCK_END(1);
	LOG_VAR(Mul);
	LIBCLOCK_PRINT(1);	
}
// ADD-BY-LEETEN 10/28/2012-END

int
main(int argn, char* argv[])
{
	size_t uNrOfElements = 100000;
	#if	0	// DEL-BY-LEETEN 10/28/2012-BEGIN
	size_t uNrOfTests = 10000;

	LIBCLOCK_INIT(1, "main");
	#endif		// DEL-BY-LEETEN 10/28/2012-END

	// test int
	// DEL-BY-LEETEN 10/28/2012:	vector<int> vi1;
	vi1.resize(uNrOfElements);
	// DEL-BY-LEETEN 10/28/2012:	vector<int> vi2;
	vi2.resize(uNrOfElements);
	for(size_t e = 0; e < uNrOfElements; e++)
	{
		#if	0	// MOD-BY-LEETEN 10/28/2012-FROM:
		vi1[e] = rand() % 128;
		vi2[e] = rand() % 128;
		#else		// MOD-BY-LEETEN 10/28/2012-TO:
		vi1[e] = 1;
		vi2[e] = 1;
		#endif		// MOD-BY-LEETEN 10/28/2012-END
	}
	#if	0	// MOD-BY-LEETEN 10/28/2012-FROM:
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

	// ADD-BY-LEETEN 10/28/2012-BEGIN
	int iMul;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		iMul = _Multiply<int>(vi1);	
	LIBCLOCK_END(1);
	LOG_VAR(iMul);
	// ADD-BY-LEETEN 10/28/2012-END

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

	// ADD-BY-LEETEN 10/28/2012-BEGIN
	unsigned int uiMul;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		uiMul = _Multiply<unsigned int>(vui1);	
	LIBCLOCK_END(1);
	LOG_VAR(uiMul);
	// ADD-BY-LEETEN 10/28/2012-END

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

	// ADD-BY-LEETEN 10/28/2012-BEGIN
	long lMul;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		lMul = _Multiply<long>(vl1);	
	LIBCLOCK_END(1);
	LOG_VAR(lMul);
	// ADD-BY-LEETEN 10/28/2012-END

	// test unsigned long
	vector<unsigned long> vul1;		_Copy<unsigned long, int>(vul1, vi1);
	vector<unsigned long> vul2;		_Copy<unsigned long, int>(vul2, vi2);
	// MOD-BY-LEETEN 10/28/2012-FROM:	long ulProd;
	unsigned long ulProd;
	// MOD-BY-LEETEN 10/28/2012-END
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		ulProd = _DotProduct<unsigned long>(vul1, vul2);	
	LIBCLOCK_END(1);
	LOG_VAR(ulProd);

	// MOD-BY-LEETEN 10/28/2012-FROM:	long ulSum;
	unsigned long ulSum;
	// MOD-BY-LEETEN 10/28/2012-END
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		ulSum = _Sum<unsigned long>(vul1, vul2);	
	LIBCLOCK_END(1);
	LOG_VAR(ulSum);

	// ADD-BY-LEETEN 10/28/2012-BEGIN
	unsigned long ulMul;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		ulMul = _Multiply<unsigned long>(vul1);	
	LIBCLOCK_END(1);
	LOG_VAR(ulMul);
	// ADD-BY-LEETEN 10/28/2012-END

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

	// ADD-BY-LEETEN 10/28/2012-BEGIN
	size_t uMul;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		uMul = _Multiply<size_t>(vu1);	
	LIBCLOCK_END(1);
	LOG_VAR(uMul);
	// ADD-BY-LEETEN 10/28/2012-END

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

	// ADD-BY-LEETEN 10/28/2012-BEGIN
	float fMul;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		fMul = _Multiply<float>(vf1);	
	LIBCLOCK_END(1);
	LOG_VAR(fMul);
	// ADD-BY-LEETEN 10/28/2012-END

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

	// ADD-BY-LEETEN 10/28/2012-BEGIN
	double dMul;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		dMul = _Multiply<double>(vd1);	
	LIBCLOCK_END(1);
	LOG_VAR(dMul);
	// ADD-BY-LEETEN 10/28/2012-END

	LIBCLOCK_PRINT(1);
	#else	// MOD-BY-LEETEN 10/28/2012-TO:
	_TestType<int>();
	_TestType<unsigned int>();
	_TestType<long>();
	_TestType<unsigned long>();
	_TestType<size_t>();
	_TestType<float>();
	_TestType<double>();
	#endif	// MOD-BY-LEETEN 10/28/2012-END

	return 0;
}
