#pragma once

#if	WITH_NETCDF
#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
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
#else	// MOD-BY-LEETEN 01/02/2013-TO:
#include "NetCDFBase.h"
#endif	// MOD-BY-LEETEN 01/02/2013-END
#endif	// #if	WITH_NETCDF

namespace WaveletSAT
{
	//! The class that contains variables to store WaveletSAT coefficients in NetCDF format
	/*
	For each coefficient, the value and ID of all bins are store in a 1D pool. The offset to the pool and the #non-zero bins are store in an D-dim array.
	*/
	#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
	class CSATSepDWTNetCDF
		:virtual public CBase	// ADD-BY-LEETEN 12/30/2012
	#else	// MOD-BY-LEETEN 01/02/2013-TO:
	class CSATSepDWTNetCDF:
		virtual public CNetCDFBase,
		virtual public CBase
	#endif	// MOD-BY-LEETEN 01/02/2013-END
	{
protected:	
		// ADD-BY-LEETEN 12/25/2012-BEGIN
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
		// ADD-BY-LEETEN 12/25/2012-END

		#if	0	// DEL-BY-LEETEN 01/02/2013-BEGIN
		//! NetCDF ID
		int iNcId;

		int iDeflateLevel; // ADD-BY-LEETEN 12/16/2012
		#endif	// DEL-BY-LEETEN 01/02/2013-END

		const char* szDimValue;
		int ncDimValue;

		#if	0	// DEL-BY-LEETEN 01/02/2013-BEGIN
		const char* szDimBin;
		int ncDimBin;

		const char* szDimDim;
		int ncDimDim;
		#endif	// DEL-BY-LEETEN 01/02/2013-END
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
			// MOD-BY-LEETEN 01/02/2013-FROM:	iDeflateLevel(0), // ADD-BY-LEETEN 12/16/2012
			CNetCDFBase(),
			// MOD-BY-LEETEN 01/02/2013-END
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
			#if	0	// DEL-BY-LEETEN 01/02/2013-BEGIN
			szDimBin("BIN"),
			szDimDim("DIM"),
			#endif	// DEL-BY-LEETEN 01/02/2013-END
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
