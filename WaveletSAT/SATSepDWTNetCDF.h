#pragma once

#if	WITH_NETCDF
#include <netcdf.h>
#include "lognc.h"
#endif	// #if	WITH_NETCDF

namespace WaveletSAT
{
	//! The class that contains variables for NetCDF
	class CSATSepDWTNetCDF
	{
protected:	
		//! NetCDF ID
		int iNcId;

		int iDeflateLevel; // ADD-BY-LEETEN 12/16/2012

		const char* szDimValue;
		int ncDimValue;

		const char* szDimBin;
		int ncDimBin;

		const char* szDimDim;
		int ncDimDim;

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
		const char* szVarHeaderOffset;
		int ncVarHeaderOffset;
		nc_type typeHeaderOffset;

		const char* szVarHeaderCount;
		int ncVarHeaderCount;
		nc_type typeHeaderCount;

		const char* szVarCoefBin;
		int ncVarCoefBin;
		nc_type typeCoefBin;

		const char* szVarCoefValue;
		int ncVarCoefValue;
		nc_type typeCoefValue;
public:
		// ADD-BY-LEETEN 12/16/2012-BEGIN
		enum EParameter
		{
			PARAMETER_BEGIN = 0x0A00,
			// ADD-BY-LEETEN 12/12/2012-BEGIN
			DEFLATE_LEVEL,
			// ADD-BY-LEETEN 12/12/2012-END
			PARAMETER_END
		};

		virtual	
		void
		_SetLong(
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
		// ADD-BY-LEETEN 12/16/2012-END

		CSATSepDWTNetCDF():
		iDeflateLevel(0), // ADD-BY-LEETEN 12/16/2012
			szDimValue("VALUE"),
			szDimBin("BIN"),
			szDimDim("DIM"),
			szVarHeaderOffset("HEADER_OFFSET"),
			szVarHeaderCount("HEADER_COUNT"),
			szVarCoefBin("COEF_BIN"),
			szVarCoefValue("COEF_VALUE")
		{
			pszDimTypes[0] = "COEF";
			pszDimTypes[1] = "DATA";
			pszDimTypes[2] = "LEVEL";
		};
	};
}
