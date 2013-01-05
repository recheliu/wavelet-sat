#pragma once

#include <vector>
#include <algorithm>	
using namespace std;
#include <math.h>

#include "DecoderBase.h"
#include "SATFileNetCDF.h"

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
	//! The class that load SAT from NetCDF files .
	#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
	template<
		typename DT	//!< Type of the data
	>
	#else	// MOD-BY-LEETEN 01/03/2013-TO:
	template<
		// typename DT,				//!< Type of the data
		typename ST = typeSum,		//!< Type of the sum
		typename BT = typeBin		//!< Type of the bin
	>
	#endif	// MOD-BY-LEETEN 01/03/2013-END
	class CSATFileDecoder:
		#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
		virtual public CHeaderBase,
		virtual public CSATFileNetCDF,
		virtual public CDecoderBase<DT>
		#else	// MOD-BY-LEETEN 01/03/2013-TO:
		// MOD-BY-LEETEN 01/04/2013-FROM:		virtual public CHeaderBase<ST, BT>,
		virtual public CHeaderBase,
		// MOD-BY-LEETEN 01/04/2013-END
		virtual public CSATFileNetCDF,
		virtual public CDecoderBase<ST, BT>
		#endif	// MOD-BY-LEETEN 01/03/2013-END
	{
protected:
	  bool bIsOutOfCore; // ADD-BY-LEETEN 01/02/2013
		double* pdSAT;
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0C00,
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
			CHeaderBase::_SetInteger(eName, lValue);
			// MOD-BY-LEETEN 01/03/2013-FROM:			CDecoderBase<DT>::_SetInteger(eName, lValue);
			CDecoderBase<ST, BT>::_SetInteger(eName, lValue);
			// MOD-BY-LEETEN 01/03/2013-END
			switch(eName)
			{
			#if	0	// DEL-BY-LEETEN 01/05/2012-BEGIN
			#if 0 // MOD-BY-LEETEN 01/02/2013-FROM:
			case RESET_IO_COUNTERS:
				uMinNrOfIORequest = uDataSize;
			#else // MOD-BY-LEETEN 01/02/2013-TO:
			case CDecoderBase<ST, BT>::RESET_IO_COUNTERS:
				this->uMinNrOfIORequest = uDataSize;
			#endif // MOD-BY-LEETEN 01/02/2013-END
				break;
			#endif	// DEL-BY-LEETEN 01/05/2012-END
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
			CHeaderBase::_GetInteger(eName, plValue);
			// MOD-BY-LEETEN 01/03/2013-FROM:			CDecoderBase<DT>::_GetInteger(eName, plValue);
			CDecoderBase<ST, BT>::_GetInteger(eName, plValue);
			// MOD-BY-LEETEN 01/03/2013-END
		}

		virtual
		void
		_Allocate(
			void *_Reserved = NULL
		) 
		{
		  // ADD-BY-LEETEN 01/02/2013-BEGIN
		  size_t uNrOfBins = this->UGetNrOfBins();
		  bIsOutOfCore = ( uNrOfBins * uDataSize * sizeof(*pdSAT) <= uSizeOfFullArrays )?false:true;
		  if( !bIsOutOfCore )
		    {
			pdSAT = new double[uNrOfBins * uDataSize];
			ASSERT_NETCDF(
				nc_get_var(iNcId, ncVarSAT, (void*)&pdSAT[0]) );
		    }
		  else
		    // ADD-BY-LEETEN 01/02/2013-END
			pdSAT = new double[uDataSize];
		}
		
		// ADD-BY-LEETEN 12/23/2012-BEGIN
		virtual	
		void
		_ShowStatistics
		(
			void *_Reserved = NULL
		)
		{
		  // ADD-BY-LEETEN 01/02/2013-BEGIN
		  const char *szFilepath = this->szFilepath;
		  // ADD-BY-LEETEN 01/02/2013-END
			#if		WITH_BOOST
			size_t uFileSize = fs::file_size( szFilepath );
			#else	// #if WITH_BOOST
			FILE *fp;
			fp = fopen(szFilepath, "rb");
			fseek(fp, 0, SEEK_END);
			size_t uFileSize = ftell(fp);
			fclose(fp);
			#endif	// #if WITH_BOOST
			LOG_VAR(uFileSize);
		}

		virtual 
		void
		_LoadFile
		(
			const char* szFilepath,
			void *_Reserved = NULL
		)
		{
			this->szFilepath = szFilepath;

			#if !WITH_NETCDF4 
			ASSERT_NETCDF(nc_open(
    				szFilepath,
    				NC_NOWRITE,
    				&iNcId));
			#else	// #if !WITH_NETCDF4
			ASSERT_NETCDF(nc_open(
    				szFilepath,
    				NC_NOWRITE | NC_NETCDF4,
    				&iNcId));
			#endif // #if !WITH_NETCDF4

			size_t uNrOfDims;
			ASSERT_NETCDF(nc_inq_dimid (
				iNcId,
				szDimDim,
				&ncDimDim));
			ASSERT_NETCDF(nc_inq_dimlen(
				iNcId,
				ncDimDim,
				&uNrOfDims));

			size_t uNrOfBins;
			ASSERT_NETCDF(nc_inq_dimid (
				iNcId,
				szDimBin,
				&ncDimBin));
			ASSERT_NETCDF(nc_inq_dimlen(
				iNcId,
				ncDimBin,
				&uNrOfBins));
			vector<size_t> vuDimLengths;
			for(size_t d = 0; d < uNrOfDims; d++)
			{
				char szDimName[NC_MAX_NAME+1];
				sprintf(szDimName, "DATA_DIM_%d", (unsigned int)d);
				int ncDim;
				ASSERT_NETCDF(nc_inq_dimid (
					iNcId,
					szDimName,
					&ncDim));
				vncDimData.push_back(ncDim);

				size_t uDimLength;
				ASSERT_NETCDF(nc_inq_dimlen(
					iNcId,
					ncDim,
					&uDimLength));

				vuDimLengths.push_back(uDimLength);
			}

			ASSERT_NETCDF(
				nc_inq_varid(iNcId, szVarSAT, &ncVarSAT) );

			_Set(vuDimLengths, uNrOfBins);
			_Allocate();
		}

		//! Return the sum of all bins at the given position
		virtual	
		void
		_GetAllSums
		(
			const vector<size_t>& vuPos,
			// MOD-BY-LEETEN 01/03/2013-FROM:			vector<DT>& vSums,
			vector<ST>& vSums,
			// MOD-BY-LEETEN 01/03/2013-END
			void *_Reserved = NULL
		)
		{
			if( UGetNrOfBins() != vSums.size() )
				vSums.resize(UGetNrOfBins());

			double* pdSums = new double[UGetNrOfBins()];

			size_t puStarts[NC_MAX_DIMS];
			size_t puCounts[NC_MAX_DIMS];
			for(size_t d = 0; d < UGetNrOfDims() + 1; d++)
			{
				if(0 == d)
				{
					puStarts[d] = 0;
					puCounts[d] = UGetNrOfBins();
				}
				else
				{
					puStarts[d] = vuPos[UGetNrOfDims() - d];
					puCounts[d] = 1;
				}
			}
			ASSERT_NETCDF(
				nc_get_vara_double(iNcId, ncVarSAT, puStarts, puCounts, pdSums) );

			for(size_t b = 0; b < UGetNrOfBins(); b++)
				// MOD-BY-LEETEN 01/03/2013-FROM:				vSums[b] = (DT)pdSums[b];
				vSums[b] = (ST)pdSums[b];
				// MOD-BY-LEETEN 01/03/2013-END

			delete [] pdSums;
		}

		// ADD-BY-LEETEN 01/05/2013-BEGIN
		//! Return the sum of all bins at the given position
		virtual	
		void
		_GetRegionSums
		(
			const vector<size_t>& vuLeft,
			const vector<size_t>& vuRight,
			vector<ST>& vdSums,
			void *_Reserved = NULL
		)
		{
			vdSums.assign(UGetNrOfBins(), (ST)0);

			size_t uNrOfQueries = (size_t)1 << this->UGetNrOfDims();
			vector<size_t> vuQueryPos;
			vuQueryPos.resize(this->UGetNrOfDims());
			vector<ST> vQuerySums;
			for(size_t q = 0; q < uNrOfQueries; q++)
			{
				int iSign = 1;
				for(size_t d = 0, j = q; d < this->UGetNrOfDims(); d++, j /= 2)
				{
					vuQueryPos[d] = (j % 2)?vuLeft[d]:vuRight[d];
					iSign *= (j % 2)?(-1):(+1);
				}
				_GetAllSums(vuQueryPos, vQuerySums);
				for(size_t b = 0; b < vdSums.size(); b++)
					vdSums[b] += (ST)iSign * vQuerySums[b];
			}
		}
		// ADD-BY-LEETEN 01/05/2013-END

		virtual
		void
		_DecodeBin
		(
			#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
			unsigned short usBin,
			valarray<DT> &vSAT,
			#else	// MOD-BY-LEETEN 01/03/2013-TO:
			const BT& usBin,
			valarray<ST> &vSAT,
			#endif	// MOD-BY-LEETEN 01/03/2013-END
			void *_Reserved = NULL
		)
		{
			if( uDataSize != vSAT.size() )
				vSAT.resize(uDataSize);

			// ADD-BY-LEETEN 01/02/2013-BEGIN
		  if( !bIsOutOfCore )
		    {
		      for(size_t d = 0; d < uDataSize; d++)
				// MOD-BY-LEETEN 01/03/2013-FROM:			vSAT[d] = (DT)pdSAT[(size_t)usBin * uDataSize + d];
				vSAT[d] = (ST)pdSAT[(size_t)usBin * uDataSize + d];
				// MOD-BY-LEETEN 01/03/2013-END
		      return;
		    }
		  // ADD-BY-LEETEN 01/02/2013-END
			size_t puStarts[NC_MAX_DIMS];
			size_t puCounts[NC_MAX_DIMS];
			for(size_t d = 0; d < UGetNrOfDims() + 1; d++)
			{
				if( !d )
				{
					puStarts[d] = (size_t)usBin;
					puCounts[d] = 1;
				}
				else
				{
					puStarts[d] = 0;
					puCounts[d] = vuDimLengths[UGetNrOfDims() - d];
				}
			}
			ASSERT_NETCDF(
				nc_get_vara(iNcId, ncVarSAT, puStarts, puCounts, (void*)&pdSAT[0]) );

			for(size_t d = 0; d < uDataSize; d++)
				// MOD-BY-LEETEN 01/03/2013-FROM:				vSAT[d] = (DT)pdSAT[d];
				vSAT[d] = (ST)pdSAT[d];
				// MOD-BY-LEETEN 01/03/2013-END
		}

		virtual
		void
		_ClampToDataSize(
			#if	0	// MOD-BY-LEETEN 01/03/2013-FROM:
			const valarray<DT>& vCoefField,
			valarray<DT>& vDataField,
			#else	// MOD-BY-LEETEN 01/03/2013-TO:
			const valarray<ST>& vCoefField,
			valarray<ST>& vDataField,
			#endif	// MOD-BY-LEETEN 01/03/2013-END
			void* _Reserved = NULL
			)
		{
			// only keep the entropy field within the data range
			if( uDataSize != vDataField.size() )
				vDataField.resize(uDataSize);
			vDataField = vCoefField;
		}

		virtual
		void
		_ClampBorder(
			// MOD-BY-LEETEN 01/03/2013-FROM:			valarray<DT>& vField,
			valarray<ST>& vField,
			// MOD-BY-LEETEN 01/03/2013-END
			const vector<int>& viLeft, 
			const vector<int>& viRight, 
			void* _Reserved = NULL
			)
		{
			vector<size_t> vuSub;
			for(size_t d = 0; d < uDataSize; d++)
			{
				vector<size_t> vuSub;
				_ConvertIndexToSub(d, vuSub, vuDimLengths);
				bool bIsNearBorder = false;
				for(size_t dim = 0; dim < this->UGetNrOfDims(); dim++)
					if( 0 > (int)vuSub[dim] + viLeft[dim] || 
							(int)vuSub[dim] + viLeft[dim] >= vuDimLengths[dim] ||
						0 > (int)vuSub[dim] + viRight[dim] || 
							(int)vuSub[dim] + viRight[dim] >= vuDimLengths[dim] )
					{
						bIsNearBorder = true;
						break;
					}

				if( bIsNearBorder )
				{
				  #if 0 // MOD-BY-LEETEN 01/04/2013-FROM:
					vField[d] = (DT)0;
					#else // MOD-BY-LEETEN 01/04/2013-TO:
					vField[d] = (ST)0;
					#endif // MOD-BY-LEETEN 01/04/2013-END
					continue;
				}
			}
		}

		virtual
		~CSATFileDecoder()
		{
			if( pdSAT )
			{
				delete [] pdSAT;
				pdSAT = NULL;
			}
		}
	};
}
