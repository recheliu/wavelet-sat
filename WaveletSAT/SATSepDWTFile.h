#pragma once

#define WITH_VECTORS_FOR_COUNTED_COEFS	0

// ADD-BY-LEETEN 10/10/2012-BEGIN
#define	WITH_BOUNDARY_AWARE_DWT		0
// ADD-BY-LEETEN 10/10/2012-END

// ADD-BY-LEETEN 11/11/2012-BEGIN
#define WITH_COEF_POOL			1
// ADD-BY-LEETEN 11/11/2012-END

#include <map>	

#include <vector>
using namespace std;
#include <math.h>

#include "utils.h"	
#include "WaveletSATEncoder.h"

#include "liblog.h"	
#include "libbuf.h"

#if	WITH_NETCDF
#include <netcdf.h>
#include "lognc.h"
#endif	// #if	WITH_NETCDF

namespace WaveletSAT
{
	template<typename T>
	// The class that load the coefficients from files (in NetCDF format)
	class CSATSepDWTFile:
#if 0 // MOD-BY-LEETEN 12/16/2012-FROM:
		public CWaveletSATEncoder<T>
#else // MOD-BY-LEETEN 12/16/2012-TO:
#ifdef 	WIN32
		public CWaveletSATEncoder<T>
#else	// #ifdef WIN32
		public CWaveletSATEncoder<double>
#endif	// #ifdef WIN32
#endif // MOD-BY-LEETEN 12/16/2012-END
	{
protected:	
public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/

		enum EParameter
		{
			PARAMETER_BEGIN = 0x0900,
			PARAMETER_END
		};

		virtual 
		void
		_LoadFile
		(
			const char* szFilepath,
			void *_Reserved = NULL
		)
		{
			#if WITH_NETCDF 

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
			#endif	// #if !WITH_NETCDF4

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

			ASSERT_NETCDF(nc_inq_dimid (
				iNcId,
				szDimValue,
				&ncDimValue));

			vector<size_t> vuDimLengths;
			for(size_t t = 0; t < NR_OF_DIM_TYPES; t++)
				for(size_t d = 0; d < uNrOfDims; d++)
				{
					char szDimName[NC_MAX_NAME+1];
					sprintf(szDimName, "%s_DIM_%d", pszDimTypes[t], (unsigned int)d);
					int ncDim;
					ASSERT_NETCDF(nc_inq_dimid (
						iNcId,
						szDimName,
						&ncDim));
					vncDims.push_back(ncDim);

					size_t uDimLength;
					ASSERT_NETCDF(nc_inq_dimlen(
						iNcId,
						ncDim,
						&uDimLength));

					switch(t)
					{
					case DIM_TYPE_DATA:
						vuDimLengths.push_back(uDimLength);
						break;
					case DIM_TYPE_COEF:
						piCoefDimIds[uNrOfDims - 1 - d] = ncDim;
						break;
					}
				}

			_Set(vuDimLengths, uNrOfBins);
			_Allocate();

			//
			// now define the variable for the coef headers
			ASSERT_NETCDF(nc_inq_varid(
				iNcId,
				szVarHeaderCount,
				&ncVarHeaderCount));
			ASSERT_NETCDF(nc_inq_varid(
				iNcId,
				szVarHeaderOffset,
				&ncVarHeaderOffset));
			ASSERT_NETCDF(nc_inq_varid(
				iNcId,
				szVarCoefValue,
				&ncVarCoefValue));
			ASSERT_NETCDF(nc_inq_varid(
				iNcId,
				szVarCoefBin,
				&ncVarCoefBin));

			/////////////////////////////////////////////////////////////////////
			// define the #non-zero bins per coef.
			#if !WITH_NETCDF4
			TBuffer<int> pHeaderOffset;
			typeHeaderOffset = NC_INT;
			TBuffer<int> pHeaderCount;
			typeHeaderCount = NC_INT;
			TBuffer<int> pCoefBin;
			typeCoefBin = NC_INT;
			#else // #if !WITH_NETCDF4
			TBuffer<unsigned long long> pHeaderOffset;
			#if	0	// MOD-BY-LEETEN 12/16/2012-FROM:
			typeHeaderOffset = NC_ULONGLONG;
			#else // MOD-BY-LEETEN 12/16/2012-TO:
			typeHeaderOffset = NC_UINT64;
			#endif // MOD-BY-LEETEN 12/16/2012-END
			TBuffer<unsigned int> pHeaderCount;
			typeHeaderCount = NC_UINT;
			TBuffer<unsigned int> pCoefBin;
			typeCoefBin = NC_UINT;
			#endif // #if !WITH_NETCDF4
			TBuffer<double> pCoefValue;
			typeCoefValue = NC_DOUBLE;

			// convert the basis id to its level subscript
			vector<size_t> vuBasisHeaderSize;
			vuBasisHeaderSize.resize(UGetNrOfDims());
			vector<size_t> vuBasisHeaderBase;
			vuBasisHeaderBase.resize(UGetNrOfDims());

			size_t puStart[NC_MAX_DIMS];
			size_t puCount[NC_MAX_DIMS];
			for(size_t c = 0; c < this->uNrOfUpdatingCoefs; c++)
			{
				vector<size_t> vuLevelSub;
				_ConvertIndexToSub(c, vuLevelSub, this->vuDimMaxLevels);

				size_t uNrOfBasisHeader = 1;
				for(size_t d = 0; d < vuLevelSub.size(); d++)
				{
					// From this subscript, we can get the dim. length in this basis. 
					size_t uLevel = vuLevelSub[d];
					size_t uLen = (!uLevel)?1:(1<<(uLevel - 1));

					vuBasisHeaderSize[d] = uLen;
					uNrOfBasisHeader *= uLen;

					// we can also decide it base in the n-dim. pool
					vuBasisHeaderBase[d] = (!uLevel)?0:uLen;

					puStart[UGetNrOfDims() - 1 - d] = vuBasisHeaderBase[d];
					puCount[UGetNrOfDims() - 1 - d] = vuBasisHeaderSize[d];
				}
				pHeaderOffset.alloc(uNrOfBasisHeader);
				pHeaderCount.alloc(uNrOfBasisHeader);

				// read the header for the current basis
				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarHeaderOffset,
					puStart,
					puCount,
					(void*)&pHeaderOffset[0]));

				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarHeaderCount,
					puStart,
					puCount,
					(void*)&pHeaderCount[0]));

				// compute the total number of coefficients in the current basis
				size_t uNrOfBasisCoefs = 0;
				for(size_t h = 0; h < uNrOfBasisHeader; h++)
					uNrOfBasisCoefs += (size_t)pHeaderCount[h];
				puStart[0] = (size_t)pHeaderOffset[0];
				puCount[0] = uNrOfBasisCoefs;

				// read the coefficients of the current basis
				pCoefBin.alloc(uNrOfBasisCoefs);
				pCoefValue.alloc(uNrOfBasisCoefs);
				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarCoefBin,
					puStart,
					puCount,
					(void*)&pCoefBin[0]));

				ASSERT_NETCDF(nc_get_vara(
					iNcId,
					ncVarCoefValue,
					puStart,
					puCount,
					(void*)&pCoefValue[0]));

				// scan through all basis
				for(size_t coefi = 0, basis = 0; basis < uNrOfBasisHeader; basis++)
				{
					vector<size_t> vuBasisSub;
					_ConvertIndexToSub(basis, vuBasisSub, vuBasisHeaderSize);

					// scan through all bin
					for(size_t bini = 0; bini < pHeaderCount[basis]; bini++, coefi++)
						this->vcCoefPools[c]._AddAt(
							pCoefBin[coefi],
							vuBasisSub,
							pCoefValue[coefi]);
				}
				this->vcCoefPools[c]._Finalize(1.0); // ADD-BY-LEETEN 12/16/2012
			}

			// close the file
			ASSERT_NETCDF(nc_close(iNcId));
			#else	// #if WITH_NETCDF 
			#endif	// #if WITH_NETCDF 

			// DEL-BY-LEETEN 12/16/2012-TO: _Finalize();
		}

		CSATSepDWTFile():
			CWaveletSATEncoder<T>()
		{
			bIsFinalizedWithoutWavelet = true;
		}
	};
}
