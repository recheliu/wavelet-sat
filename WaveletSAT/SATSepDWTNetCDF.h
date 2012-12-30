#pragma once

#if	WITH_NETCDF
#include <netcdf.h>
#include "lognc.h"

// ADD-BY-LEETEN 01/22/2012-BEGIN
#if defined(WIN32)
#if	defined(_DEBUG)
	#pragma comment (lib, "zlibd.lib")
	#pragma comment (lib, "hdf5d.lib")
	#pragma comment (lib, "hdf5_hld.lib")
#else	// #if	defined(_DEBUG)
	#pragma comment (lib, "zlib.lib")
	#pragma comment (lib, "hdf5.lib")
	#pragma comment (lib, "hdf5_hl.lib")
#endif	// #if	defined(_DEBUG)
#pragma comment (lib, "libsrc.lib")
#pragma comment (lib, "libsrc4.lib")
#pragma comment (lib, "libdispatch.lib")
#endif	// #if defined(WIN32)
// ADD-BY-LEETEN 01/22/2012-END

#endif	// #if	WITH_NETCDF

namespace WaveletSAT
{
	//! The class that contains variables for NetCDF
	/*
	For each coefficient, the value and ID of all bins are store in a 1D pool. The offset to the pool and the #non-zero bins are store in an D-dim array.
	*/
	class CSATSepDWTNetCDF
		:virtual public CBase	// ADD-BY-LEETEN 12/30/2012
	{
protected:	
		// ADD-BY-LEETEN 12/25/2012-BEGIN
		#if !WITH_NETCDF4
		#if	0	// MOD-BY-LEETEN 12/29/2012-FROM:
		typedef int	TYPE_HEADER_COUNT;
		typedef int	TYPE_HEADER_OFFSET;
		#else	// MOD-BY-LEETEN 12/29/2012-TO:
		typedef int	TYPE_COEF_COUNT;
		typedef int	TYPE_COEF_OFFSET;
		#endif	// MOD-BY-LEETEN 12/29/2012-END
		typedef int	TYPE_COEF_BIN;
		#else // #if !WITH_NETCDF4
		typedef unsigned int 		TYPE_COEF_COUNT;
		typedef unsigned long long	TYPE_COEF_OFFSET;
		typedef unsigned int		TYPE_COEF_BIN;
		#endif // #if !WITH_NETCDF4
		typedef double TYPE_COEF_VALUE;
		// ADD-BY-LEETEN 12/25/2012-END

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
		#if	0	// MOD-BY-LEETEN 12/29/2012-FROM:
		const char* szVarHeaderOffset;
		int ncVarHeaderOffset;
		nc_type typeHeaderOffset;

		const char* szVarHeaderCount;
		int ncVarHeaderCount;
		nc_type typeHeaderCount;
		#else	// MOD-BY-LEETEN 12/29/2012-TO:
		const char* szVarCoefOffset;
		int ncVarCoefOffset;
		nc_type typeCoefOffset;

		const char* szVarCoefCount;
		int ncVarCoefCount;
		nc_type typeCoefCount;
		#endif	// MOD-BY-LEETEN 12/29/2012-END
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
		// ADD-BY-LEETEN 12/16/2012-END

		CSATSepDWTNetCDF():
		iDeflateLevel(0), // ADD-BY-LEETEN 12/16/2012
			// ADD-BY-LEETEN 12/25/2012-BEGIN
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
			// ADD-BY-LEETEN 12/25/2012-END
			szDimValue("VALUE"),
			szDimBin("BIN"),
			szDimDim("DIM"),
			#if	0	// MOD-BY-LEETEN 12/29/2012-FROM:
			szVarCoefOffset("HEADER_OFFSET"),
			szVarCoefCount("HEADER_COUNT"),
			#else	// MOD-BY-LEETEN 12/29/2012-TO:
			szVarCoefOffset("COEF_OFFSET"),
			szVarCoefCount("COEF_COUNT"),
			#endif	// MOD-BY-LEETEN 12/29/2012-END
			szVarCoefBin("COEF_BIN"),
			szVarCoefValue("COEF_VALUE")
		{
			pszDimTypes[0] = "COEF";
			pszDimTypes[1] = "DATA";
			pszDimTypes[2] = "LEVEL";
		};
	};
}
