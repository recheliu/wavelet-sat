#pragma once

#if	WITH_NETCDF
#include "NetCDFBase.h"
#endif	// #if	WITH_NETCDF

namespace WaveletSAT
{
	//! The class that contains variables to store WaveletSAT coefficients in NetCDF format
	/*
	For each coefficient, the value and ID of all bins are store in a 1D pool. The offset to the pool and the #non-zero bins are store in an D-dim array.
	*/
	class CSATSepDWTNetCDF:
		virtual public CNetCDFBase,
		virtual public CBase
	{
protected:	
		const char* szDimValue;
		int ncDimValue;

		enum {
			DIM_TYPE_COEF,
			DIM_TYPE_DATA, 
			DIM_TYPE_LEVEL, 
			NR_OF_DIM_TYPES
		};
		const char *pszDimTypes[NR_OF_DIM_TYPES];
		vector<int> vncDims;

		int piCoefDimIds[NC_MAX_DIMS];

		// Variable IDs
		const char* szVarCoefOffset;
		int ncVarCoefOffset;
		nc_type typeCoefOffset;

		const char* szVarCoefCount;
		int ncVarCoefCount;
		nc_type typeCoefCount;

		const char* szVarCoefBin;
		int ncVarCoefBin;
		nc_type typeCoefBin;

		const char* szVarCoefValue;
		int ncVarCoefValue;
		nc_type typeCoefValue;
public:
		enum EParameter
		{
			PARAMETER_BEGIN = 0x0A00,
			DEFLATE_LEVEL,
			PARAMETER_END
		};

		#if !WITH_NETCDF4
		typedef int	TYPE_COEF_COUNT;
		typedef int	TYPE_COEF_OFFSET;
		typedef int	TYPE_COEF_BIN;
		#else // #if !WITH_NETCDF4
		typedef unsigned int 		TYPE_COEF_COUNT;
		typedef unsigned long long	TYPE_COEF_OFFSET;
		typedef unsigned int		TYPE_COEF_BIN;
		#endif // #if !WITH_NETCDF4
		typedef double TYPE_COEF_VALUE;

		virtual	
		void
		_SetInteger(
			int eName,
			long lValue,
			void* _Reserved = NULL
		)
		{
		  switch(eName)
		    {
		    case DEFLATE_LEVEL:
		      iDeflateLevel = lValue;
		      break;
		    }
		}

		CSATSepDWTNetCDF():
			CNetCDFBase(),
			#if !WITH_NETCDF4
			typeCoefOffset(NC_INT),
			typeCoefCount(NC_INT),
			typeCoefBin(NC_INT),
			#else // #if !WITH_NETCDF4
			typeCoefOffset(NC_UINT64),
			typeCoefCount(NC_UINT),
			typeCoefBin(NC_UINT),
			#endif // #if !WITH_NETCDF4
			typeCoefValue(NC_DOUBLE),
			szDimValue("VALUE"),
			szVarCoefOffset("COEF_OFFSET"),
			szVarCoefCount("COEF_COUNT"),
			szVarCoefBin("COEF_BIN"),
			szVarCoefValue("COEF_VALUE")
		{
			pszDimTypes[0] = "COEF";
			pszDimTypes[1] = "DATA";
			pszDimTypes[2] = "LEVEL";
		};
	};
}

/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
