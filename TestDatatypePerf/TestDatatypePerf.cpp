#include <typeinfo>	// ADD-BY-LEETEN 10/28/2012
#include<vector>
using namespace std;
#include <stdio.h>
#include <assert.h>
#include <stdlib.h> // ADD-BY-LEETEN 10/22/2012
#include <string.h> // ADD-BY-LEETEN 10/29/2012 
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

	// ADD-BY-LEETEN 01/22/2012-BEGIN
	// Allocate buffer to hold the result. This forces the compiler to compute all iterations
	vector<T> vProds;
	vProds.resize(uNrOfTests);
	vector<T> vSums;
	vSums.resize(uNrOfTests);
	vector<T> vMuls;
	vMuls.resize(uNrOfTests);
	// ADD-BY-LEETEN 01/22/2012-END

	LIBCLOCK_END(1);

	LIBCLOCK_BEGIN(1);	
	T Prod;
	for(size_t t = 0; t < uNrOfTests; t++)
		vProds[t] = _DotProduct<T>(v1, v2);	
	Prod = vProds[0];
	LIBCLOCK_END(1);
	LOG_VAR(Prod);

	T Sum;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		vSums[t] = _Sum<T>(v1, v2);	
	Sum = vSums[0];
	LIBCLOCK_END(1);
	LOG_VAR(Sum);

	T Mul;
	LIBCLOCK_BEGIN(1);	
	for(size_t t = 0; t < uNrOfTests; t++)
		vMuls[t] = _Multiply<T>(v1);	
	Mul = vMuls[0];
	LIBCLOCK_END(1);
	LOG_VAR(Mul);
	LIBCLOCK_PRINT(1);	
}
// ADD-BY-LEETEN 10/28/2012-END

int
main(int argn, char* argv[])
{
	size_t uNrOfElements = 100000;
	// test int
	vi1.resize(uNrOfElements);
	vi2.resize(uNrOfElements);
	for(size_t e = 0; e < uNrOfElements; e++)
	{
		vi1[e] = 1;
		vi2[e] = 1;
	}
	_TestType<int>();
	_TestType<unsigned int>();
	_TestType<long>();
	_TestType<unsigned long>();
	_TestType<size_t>();
	_TestType<float>();
	_TestType<double>();

	return 0;
}
