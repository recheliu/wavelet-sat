#pragma once

#if	WITH_NETCDF
#include "NetCDFBase.h"
#endif	// #if	WITH_NETCDF

namespace WaveletSAT
{
	//! The class that contains variables to store the SAT in NetCDF format
	/*
	For each coefficient, the value and ID of all bins are store in a 1D pool. The offset to the pool and the #non-zero bins are store in an D-dim array.
	*/
	class CSATFileNetCDF:
		virtual public CNetCDFBase,
		virtual public CBase	
	{
protected:	
		vector<int> vncDimData;

		int ncVarSAT;
		const char* szVarSAT;

		const char* szVarCoefOffset;
		int			ncVarCoefOffset;
		nc_type		nctypeCoefOffset;

		const char* szVarCoefCount;
		int			ncVarCoefCount;
		nc_type		nctypeCoefCount;

		int			ncDimCoef;
		const char* szDimCoef;

		const char* szVarCoefBin;
		int			ncVarCoefBin;
		nc_type		nctypeCoefBin;

		const char* szVarCoefSum;
		int			ncVarCoefSum;
		nc_type		nctypeCoefSum;

public:
		static const char* SZGetVarCoefOffset()	{	return	"COEF_OFFSET";	}
		static const char* SZGetVarCoefCount()	{	return	"COEF_COUNT";	}
		static const char* SZGetVarCoefBin()	{	return	"COEF_BIN";		}
		static const char* SZGetVarCoefSum()	{	return	"COEF_SUM";		}

		typedef		unsigned long long typeCoefOffset;
		typedef		unsigned short typeCoefCount;

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0D00,
			PARAMETER_END
		};

		virtual	
		void
		_SetInteger(
			int eName,
			long lValue,
			void* _Reserved = NULL
		)
		{
			CNetCDFBase::_SetInteger(eName, lValue);
		}

		CSATFileNetCDF():
			szVarSAT("SAT"),
			szDimCoef("COEF"),
			szVarCoefOffset(SZGetVarCoefOffset()),
			szVarCoefCount(	SZGetVarCoefCount()),
			szVarCoefBin(	SZGetVarCoefBin()),
			szVarCoefSum(	SZGetVarCoefSum()),
			nctypeCoefOffset(	NC_UINT64),
			nctypeCoefCount(	NC_USHORT),
			nctypeCoefBin(		NC_USHORT),
			nctypeCoefSum(		NC_DOUBLE),
			CNetCDFBase()
		{
		};
	};
}

/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
