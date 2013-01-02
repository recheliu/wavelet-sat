#pragma once

#include <vector>
#include <algorithm>	
using namespace std;
#include <math.h>

#include "Base.h"

#include "liblog.h"	
#if	!WITH_SMART_PTR	// ADD-BY-LEETEN 12/30/2012
#include "libbuf.h"
// ADD-BY-LEETEN 12/30/2012-BEGIN
#else	// #if	!WITH_SMART_PTR
#include <boost/shared_array.hpp>
#endif	// #if	!WITH_SMART_PTR
// ADD-BY-LEETEN 12/30/2012-END

namespace WaveletSAT
{
	template<
		typename DT	//!< Type of the data
	>
	class CDecoderBase:
		virtual public CBase
	{
	protected:
			//! The accumulated #I/O requests since the I/O counters are reset.
			size_t uAccumNrOfIORequest;

			//! The max #I/O requests per query since the I/O counters are reset.
			size_t uMaxNrOfIORequest;

			//! The min #I/O requests per query since the I/O counters are reset.
			size_t uMinNrOfIORequest;

			//! The #query since the I/O counters are reset.
			size_t uNrOfQueries;
			// ADD-BY-LEETEN 12/28/2012-END

			const char* szFilepath;	

			bool bIsPrintingDecodeBinTiming;

	public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0B00,

			ACCUM_NR_OF_IO_REQUESTS,
			MAX_NR_OF_IO_REQUESTS,
			MIN_NR_OF_IO_REQUESTS,
			RESET_IO_COUNTERS,
			
			//! Parameter specifying whether to print the timing for the method _Decodebin().
			PRINT_DECODE_BIN_TIMING,
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
			case RESET_IO_COUNTERS:
				uAccumNrOfIORequest = 0;
				uMaxNrOfIORequest = 0;
				uMinNrOfIORequest = 0;
				uNrOfQueries = 0;
				break;

			case PRINT_DECODE_BIN_TIMING:
				bIsPrintingDecodeBinTiming = (!lValue)?false:true;
				break;
			}
		}

		virtual	
		void
		_GetInteger(
			int eName,
			long *plValue,
			void* _Reserved = NULL
		)
		{
			switch(eName)
			{
			case ACCUM_NR_OF_IO_REQUESTS:
				*plValue = (long)uAccumNrOfIORequest;
				break;
			case MAX_NR_OF_IO_REQUESTS:
				*plValue = (long)uMaxNrOfIORequest;
				break;
			case MIN_NR_OF_IO_REQUESTS:
				*plValue = (long)uMinNrOfIORequest;
				break;
			}
		}

		virtual
		void
		_Allocate(
			void *_Reserved = NULL
		) = 0;
		
		// ADD-BY-LEETEN 12/23/2012-BEGIN
		virtual	
		void
		_ShowStatistics
		(
			void *_Reserved = NULL
		) = 0;

		virtual 
		void
		_LoadFile
		(
			const char* szFilepath,
			void *_Reserved = NULL
		) = 0;

		//! Return the sum of all bins at the given position
		virtual	
		void
		_GetAllSums
		(
			const vector<size_t>& vuPos,
			vector<DT>& vdSums,
			void *_Reserved = NULL
		) = 0;

		// ADD-BY-LEETEN 12/29/2012-BEGIN
		virtual
		void
		_DecodeBin
		(
			unsigned short usBin,
			valarray<DT> &vSAT,
			void *_Reserved = NULL
		) = 0;

		virtual
		void
		_ClampToDataSize(
			const valarray<DT>& vCoefField,
			valarray<DT>& vDataField,
			void* _Reserved = NULL
			) = 0;

		virtual
		void
		_ClampBorder(
			valarray<DT>& vField,
			const vector<int>& viLeft, 
			const vector<int>& viRight, 
			void* _Reserved = NULL
			) = 0;

		CDecoderBase():
			bIsPrintingDecodeBinTiming(false)
		{
		}
	};
}
