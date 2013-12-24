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

#include "Base.h"

namespace WaveletSAT
{
	//! The class that contains variables for NetCDF
	/*
	For each coefficient, the value and ID of all bins are store in a 1D pool. The offset to the pool and the #non-zero bins are store in an D-dim array.
	*/
	class CNetCDFBase:
		virtual public CBase	
	{
protected:	
		int iNcId;

		int iDeflateLevel; // ADD-BY-LEETEN 12/16/2012

		char szNetCDFFilepath[NC_MAX_NAME];

		const char* szDimBin;
		int ncDimBin;

		const char* szDimDim;
		int ncDimDim;
public:
		// ADD-BY-LEETEN 12/16/2012-BEGIN
		enum EParameter
		{
			PARAMETER_BEGIN = 0x0E00,
			DEFLATE_LEVEL,
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

		CNetCDFBase():
			iDeflateLevel(0),
			szDimBin("BIN"),
			szDimDim("DIM"),
			iNcId(0)
		{
		};
	};
}

/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
