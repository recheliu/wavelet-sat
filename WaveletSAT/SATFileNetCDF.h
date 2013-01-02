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
public:
		// ADD-BY-LEETEN 12/16/2012-BEGIN
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
		// ADD-BY-LEETEN 12/16/2012-END

		CSATFileNetCDF():
			szVarSAT("SAT"),
			CNetCDFBase()
		{
		};
	};
}
