#pragma once

#include <iostream>
#include <vector>
#include <map>
using namespace std;
#include <math.h>

#include <netcdf.h>

/*
#pragma comment (lib, "libszip.lib")
#pragma comment (lib, "zlib.lib")
#pragma comment (lib, "hdf5.lib")
#pragma comment (lib, "hdf5_cpp.lib")
#pragma comment (lib, "hdf5_hl.lib")
#pragma comment (lib, "hdf5_hl_cpp.lib")
#pragma comment (lib, "netcdf.lib")
*/

#include "liblog.h"
#include "lognc.h"

#include "EncoderBase.h"

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

	//! The class that create SAT in core.
	/*!
	*/
	template<typename DT, typename ST>
	class CSATFileEncoder:
		public CHeaderBase,
		public CEncoderBase<DT, ST>
	{
protected:	
		//! The storage to store the original data.
		vector< map<size_t, ST> > vmapHists;

		char szNetCdfPathFilename[NC_MAX_NAME];

		int iNcId;

		int piNcDimIds[NC_MAX_DIMS];

		int iNcVarId;
		
		char *pszNcDimNames[NC_MAX_DIMS];

		////////////////////////////////////////////////////////////////////
		/*
		The protected interface. 
		*/
		//! Update the specified bin.
		virtual
		void 
		_UpdateBin
		(
			const vector<size_t>& vuPos, 
			const DT& value,
			size_t uBin, 
			const ST& weight,
			void *_Reserved = NULL
		)
		{
			size_t uIndex = UConvertSubToIndex(vuPos, vuDimLengths);
			map<size_t, ST>& mapHist = vmapHists[uIndex];
			map<size_t, double>::iterator imapHist = mapHist.find(uBin);
			if(mapHist.end() == imapHist )
				mapHist.insert(pair<size_t, ST>(uBin, weight));
			else
				imapHist->second += weight;
		}

public:
		////////////////////////////////////////////////////////////////////
		/*
		The public interface. 
		*/
		enum EParameter
		{
			PARAMETER_BEGIN = 0x0300,
			PARAMETER_END
		};

		//! Finalize the computation of SAT
		/*!
		*/
		virtual 
		void 
		_Finalize
		(
			void *_Reserved = NULL
		)
		{
			/////////////////////////////////////////////////////
			// create the NetCDF file
			sprintf(szNetCdfPathFilename, "sat.nc");
			/*
			sprintf(szNetCdfPathFilename, "%s/%s.rank_%d.nc", szPath, szFilenamePrefix, iRank);
			*/

			// Create the file.
			ASSERT_NETCDF(nc_create(
    				szNetCdfPathFilename,
    				NC_CLOBBER,
    				&iNcId));

			// dimension
			int piDimLengths[NC_MAX_DIMS];
			static char* pszDefaultDimNames[] = {"X", "Y", "Z"};
			static size_t uNrOfDefaultDimNames = sizeof(pszDefaultDimNames)/sizeof(pszDefaultDimNames[0]);

			for(int d = 0; d < UGetNrOfDims() + 1; d++)
			{
				size_t uDimLength;
				char* szNcDimName;
				if( 0 == d )
				{
					uDimLength = UGetNrOfBins();
					szNcDimName = strdup("Bin");
				}
				else
				{
					size_t uDimDiff = UGetNrOfDims() - d;
					uDimLength = (int)this->vuDimLengths[uDimDiff];
					if( d < uNrOfDefaultDimNames )
						szNcDimName = strdup(pszDefaultDimNames[uDimDiff]);
					else
					{
						szNcDimName = (char*)calloc(8, 1);
						// MOD-BY-LEETEN 10/28/2012-FROM:						sprintf(szNcDimName, "DIM%d", uDimDiff);
						sprintf(szNcDimName, "DIM%d", (unsigned int)uDimDiff);
						// MOD-BY-LEETEN 10/28/2012-END
					}
				}

				piDimLengths[d] = (int)uDimLength;
				pszNcDimNames[d] = szNcDimName;
				ASSERT_NETCDF(nc_def_dim(
    							iNcId,
    							pszNcDimNames[d],
    							piDimLengths[d],
    							&piNcDimIds[d]));
			}

			// now define the variable for the bin SAT
			ASSERT_NETCDF(nc_def_var(
				   iNcId,
				   "SAT", 
				   NC_DOUBLE,
				   UGetNrOfDims() + 1,
				   piNcDimIds,
				   &iNcVarId));

			// finish the definition mode
			ASSERT_NETCDF(nc_enddef(iNcId));

			/////////////////////////////////////////
			// start to create the SAT
			vector<size_t> vuScanLineBase;	// the base coordinate of the scane line
			vuScanLineBase.resize(this->UGetNrOfDims());

			vector<size_t> vuOtherDimLengths;	// the base coordinate of the scane line
			vuOtherDimLengths.resize(this->UGetNrOfDims());

			ST *pSAT;
			pSAT = new ST[this->uDataSize];
			for(size_t 	b = 0;
					b < UGetNrOfBins(); 
					b++)
			{
				for(size_t i = 0; i < this->uDataSize; i++)
				{
					const map<size_t, ST>& mapHist = vmapHists[i];
					// MOD-BY-LEETEN 10/28/2012-FROM:					map<size_t, ST>::const_iterator imapHist = mapHist.find(b);
					typename map<size_t, ST>::const_iterator imapHist = mapHist.find(b);
					// MOD-BY-LEETEN 10/28/2012-END
					pSAT[i] = ( mapHist.end() == imapHist )?0:imapHist->second;
				}
				
				for(size_t uOffset = 1, d = 0; d < this->UGetNrOfDims(); uOffset *= vuDimLengths[d], d++)
				{
					size_t uNrOfScanLines = this->uDataSize / this->vuDimLengths[d];
					vuOtherDimLengths = vuDimLengths;
					vuOtherDimLengths[d] = 1;

					for(size_t i = 0; i < uNrOfScanLines; i++)
					{
						_ConvertIndexToSub(i, vuScanLineBase, vuOtherDimLengths);
						size_t uScanLineBase = UConvertSubToIndex(vuScanLineBase, vuDimLengths);
						for(size_t 	j = 1, uIndex = uScanLineBase; 
								j < vuDimLengths[d]; 
								j++, uIndex += uOffset)
							pSAT[uIndex + uOffset] += pSAT[uIndex];
					}
				}

				// dump this SAT
				size_t puStart[NC_MAX_DIMS];
				size_t puCount[NC_MAX_DIMS];
				for(size_t d = 0; d < UGetNrOfDims() + 1; d++)
					if( 0 == d )
					{
						puStart[d] = b;
						puCount[d] = 1;
					}
					else
					{
						puStart[d] = 0;
						puCount[d] = this->vuDimLengths[UGetNrOfDims() - d];
					}

				ASSERT_NETCDF(nc_put_vara_double(
						   iNcId,
						   iNcVarId,
						   puStart,
						   puCount,
						   &pSAT[0]));
			}
			ASSERT_NETCDF(nc_close(iNcId));
		
			delete [] pSAT;

			/////////////////////////////////////////////
			// now reopen the file in read-only mode
			ASSERT_NETCDF(nc_open(
    				szNetCdfPathFilename,
    				NC_NOWRITE,
    				&iNcId));

			ASSERT_NETCDF(
				nc_inq_varid(iNcId, "SAT", &iNcVarId) );
		}

		//! Allocate the space to store coefficients for all bins. 
		/*! 
		*/
		virtual 
		void 
		_Allocate
		(
			void *_Reserved = NULL
		)
		{
			if( uNrOfBins * uDataSize > uSizeOfFullArrays )
			{
				LOG_ERROR(cerr<<"Exceed the specified application-side memory capacity.");
				exit(EXIT_FAILURE);
			}
			vmapHists.resize(uDataSize);
		}
		
		//! Compute and display statistics for the computed wavelet coefficients.
		virtual 
		void
		_ShowStatistics
		(
			void *_Reserved = NULL
		)
		{
		}

		//! Return the sum of all bins at the given position
		virtual
		void
		_GetAllSums
		(
			const vector<size_t>& vuPos,
			vector<ST>& vdSums,
			void *_Reserved = NULL
		)
		{
			vdSums.resize(UGetNrOfBins());

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
				nc_get_vara_double(iNcId, iNcVarId, puStarts, puCounts, pdSums) );

			for(size_t b = 0; b < UGetNrOfBins(); b++)
				vdSums[b] = pdSums[b];

			delete [] pdSums;
		}		
		
		CSATFileEncoder()
		{
		}

		~CSATFileEncoder()
		{
			ASSERT_NETCDF(nc_close(iNcId));

			for(size_t d = 0; d < UGetNrOfDims() + 1; d++)
				free(pszNcDimNames[d]);
		}
	};
}
