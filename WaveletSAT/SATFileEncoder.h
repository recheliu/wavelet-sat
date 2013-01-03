#pragma once

// ADD-BY-LEETEN 12/28/2012-BEGIN
#if	WITH_BOOST
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;
#endif	// #if	WITH_BOOST
// ADD-BY-LEETEN 12/28/2012-END

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
#include "SATFileNetCDF.h"	// ADD-BY-LEETEN 01/02/2013

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	/* usage: 
	Setup #bins (to allocate #coefficients), dimensions (so the program can pre-compute the SAT of the wavelet basis), 
	*/

	//! The class that creates SAT in NetCDF files.
	/*!
	*/
	template<typename DT, typename ST>
	class CSATFileEncoder:
		virtual public CHeaderBase,
		virtual public CSATFileNetCDF,	// ADD-BY-LEETEN 01/02/2013
		virtual public CEncoderBase<DT, ST>
	{
protected:	
		//! The storage to store the original data.
		vector< map<size_t, ST> > vmapHists;

		char szNetCdfPathFilename[NC_MAX_NAME];

		#if	0	// DEL-BY-LEETEN 01/02/2013-BEGIN
		int iNcId;

		int piNcDimIds[NC_MAX_DIMS];

		int iNcVarId;
		
		char *pszNcDimNames[NC_MAX_DIMS];

		int iDeflateLevel; // ADD-BY-LEETEN 11/09/2012
		#endif	// DEL-BY-LEETEN 01/02/2013-END

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
			PARAMETER_BEGIN = 0x0600,
			
			DEFLATE_LEVEL,

			PARAMETER_END
		};

		// ADD-BY-LEETEN 11/09/2012-BEGIN
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
				iDeflateLevel = (int)lValue;
				break;
			}
			CHeaderBase::_SetInteger(eName, lValue);
			// CEncoderBase<DT, ST>::_SetInteger(eName, lValue);
		}
		// ADD-BY-LEETEN 11/09/2012-END

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
			#if !WITH_NETCDF4 
			ASSERT_NETCDF(nc_create(
    				szNetCdfPathFilename,
    				NC_CLOBBER,
    				&iNcId));
			#else	// #if !WITH_NETCDF4 
			ASSERT_NETCDF(nc_create(
    				szNetCdfPathFilename,
				NC_CLOBBER | NC_NETCDF4,
    				&iNcId));
			#endif	// #if !WITH_NETCDF4 

			#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
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
						sprintf(szNcDimName, "DIM%d", (unsigned int)uDimDiff);
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

			// ADD-BY-LEETEN 11/09/2012-BEGIN
#if WITH_NETCDF4 
			ASSERT_NETCDF(nc_def_var_deflate(
				   iNcId,
				   iNcVarId, 
				   0, 
				   1, 
				   iDeflateLevel));
#endif // #if WITH_NETCDF4
			// ADD-BY-LEETEN 11/09/2012-END
			#else	// MOD-BY-LEETEN 01/02/2013-TO:
			// define the dimension for #bins
			ASSERT_NETCDF(nc_def_dim(
						iNcId,
						szDimBin,
						(int)this->UGetNrOfBins(),
						&ncDimBin) );
			// define the dimension for #bins
			ASSERT_NETCDF(nc_def_dim(
						iNcId,
						szDimDim,
						(int)this->UGetNrOfDims(),
						&ncDimDim ) );

			size_t uNrOfDims = this->UGetNrOfDims();
			for(int d = 0; d < uNrOfDims; d++)
			{
				size_t uDimLength = (size_t)this->vuDimLengths[d];
				char szNcDimName[NC_MAX_NAME+1];
				sprintf(szNcDimName, "DATA_DIM_%d", (unsigned int)d);
				int iNcDimId;
				ASSERT_NETCDF(nc_def_dim(
    							iNcId,
    							szNcDimName,
    							(int)uDimLength,
    							&iNcDimId));
				vncDimData.push_back(iNcDimId);
			}

			// now define the variable for the bin SAT
			int piNcDimIds[NC_MAX_DIMS];
			piNcDimIds[0] = ncDimBin;
			for(size_t d = 0; d < uNrOfDims; d++)
				piNcDimIds[d + 1] = vncDimData[uNrOfDims - 1 - d];
			ASSERT_NETCDF(nc_def_var(
				   iNcId,
				   szVarSAT, 
				   NC_DOUBLE,
				   uNrOfDims + 1,
				   piNcDimIds,
				   &ncVarSAT));

			#if WITH_NETCDF4 
			ASSERT_NETCDF(nc_def_var_deflate(
				   iNcId,
				   ncVarSAT, 
				   0, 
				   1, 
				   iDeflateLevel));
			#endif // #if WITH_NETCDF4
			#endif	// MOD-BY-LEETEN 01/02/2013-END
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
					typename map<size_t, ST>::const_iterator imapHist = mapHist.find(b);
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

				#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
				ASSERT_NETCDF(nc_put_vara_double(
						   iNcId,
						   iNcVarId,
						   puStart,
						   puCount,
						   &pSAT[0]));
				#else	// MOD-BY-LEETEN 01/02/2013-TO:
				#if 0 				// MOD-BY-LEETEN 01/02/2013-FROM:
				ASSERT_NETCDF(nc_put_vara_double(
						   iNcId,
						   ncVarSAT,
						   puStart,
						   puCount,
						   &pSAT[0]));
				#else // MOD-BY-LEETEN 01/02/2013-TO:
				ASSERT_NETCDF(nc_put_vara(
						   iNcId,
						   ncVarSAT,
						   puStart,
						   puCount,
						   (void*)&pSAT[0]));
				#endif // MOD-BY-LEETEN 01/02/2013-END
				#endif	// MOD-BY-LEETEN 01/02/2013-END
			}
			ASSERT_NETCDF(nc_close(iNcId));
		
			delete [] pSAT;

			#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
			/////////////////////////////////////////////
			// now reopen the file in read-only mode
			#if !WITH_NETCDF4 
			ASSERT_NETCDF(nc_open(
    				szNetCdfPathFilename,
    				NC_NOWRITE,
    				&iNcId));
			#else	// #if !WITH_NETCDF4
			ASSERT_NETCDF(nc_open(
    				szNetCdfPathFilename,
    				NC_NOWRITE | NC_NETCDF4,
    				&iNcId));
			#endif // #if !WITH_NETCDF4
			#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
			ASSERT_NETCDF(
				nc_inq_varid(iNcId, "SAT", &iNcVarId) );
			#else	// MOD-BY-LEETEN 01/02/2013-TO:
			ASSERT_NETCDF(
				nc_inq_varid(iNcId, szVarSAT, &ncVarSAT) );
			#endif	// MOD-BY-LEETEN 01/02/2013-END
			#else	// MOD-BY-LEETEN 01/02/2013-TO:
			iNcId = 0;
			#endif	// MOD-BY-LEETEN 01/02/2013-END
		}

		// ADD-BY-LEETEN 12/12/2012-BEGIN
		//! Save the coefficients to a file
		virtual 
		void
		_SaveFile
		(
			const char* szFilepathPrefix,
			void *_Reserved = NULL
		)
		{
			#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
			#ifdef		WIN32
			const char *szCopyCommand = "copy";
			#else	//	#ifdef WIN32
			const char *szCopyCommand = "cp";
			#endif	//	#ifdef WIN32

			#if !WITH_NETCDF4
			const char* szFileExt = "nc";
			#else // #if !WITH_NETCDF4
			const char* szFileExt = "nc4";
			#endif // #if !WITH_NETCDF4

			char szCommandLine[1024];	// ADD-BY-LEETEN 01/02/2013
			sprintf(szCommandLine, "%s sat.nc %s.%s", szCopyCommand, szFilepathPrefix, szFileExt);
			// DEL-BY-LEETEN 01/02/2013:	#endif
			system(szCommandLine);
			#else	// MOD-BY-LEETEN 01/02/2013-TO:
			const char* szExt;
			#if WITH_NETCDF4
			szExt = "nc4";
			#else	// #if WITH_NETCDF4
			szExt = "nc";
			#endif	// #if WITH_NETCDF4
			char szFilepath[NC_MAX_NAME];
			sprintf(szFilepath, "%s.%s", szFilepathPrefix, szExt);
			remove(szFilepath);
			rename("sat.nc", szFilepath);
			#endif	// MOD-BY-LEETEN 01/02/2013-END
		}
		// ADD-BY-LEETEN 12/12/2012-END

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
			#if	0	// DEL-BY-LEETEN 01/02/2013-BEGIN
			if( uNrOfBins * uDataSize > uSizeOfFullArrays )
			{
				LOG_ERROR(cerr<<"Exceed the specified application-side memory capacity.");
				exit(EXIT_FAILURE);
			}
			#endif	// DEL-BY-LEETEN 01/02/2013-END
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
			// ADD-BY-LEETEN 12/28/2012-BEGIN
			#if		WITH_BOOST
			#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
			fs::path pathNativePath( szFilepath, fs::native );
			size_t uFileSize = fs::file_size( pathNativePath );
			#else	// MOD-BY-LEETEN 01/02/2013-TO:
			size_t uFileSize = fs::file_size( szNetCdfPathFilename );
			#endif	// MOD-BY-LEETEN 01/02/2013-END
			#else	// #if WITH_BOOST
			// ADD-BY-LEETEN 12/28/2012-END
		  // ADD-BY-LEETEN 11/09/2012-BEGIN
		  FILE *fp;
		  fp = fopen(szNetCdfPathFilename, "rb");
		  fseek(fp, 0, SEEK_END);
		  size_t uFileSize = ftell(fp);
		  fclose(fp);
			#endif	// #if WITH_BOOST	// ADD-BY-LEETEN 12/28/2012
		  LOG_VAR(uFileSize);
		  // ADD-BY-LEETEN 11/09/2012-END
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
			// ADD-BY-LEETEN 01/02/2013-BEGIN
			if( !iNcId )
			{
				/////////////////////////////////////////////
				// now reopen the file in read-only mode
				#if !WITH_NETCDF4 
				ASSERT_NETCDF(nc_open(
    					szNetCdfPathFilename,
    					NC_NOWRITE,
    					&iNcId));
				#else	// #if !WITH_NETCDF4
				ASSERT_NETCDF(nc_open(
    					szNetCdfPathFilename,
    					NC_NOWRITE | NC_NETCDF4,
    					&iNcId));
				#endif // #if !WITH_NETCDF4
				ASSERT_NETCDF(
					nc_inq_varid(iNcId, szVarSAT, &ncVarSAT) );
			}
			// ADD-BY-LEETEN 01/02/2013-END

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
			#if	0	// MOD-BY-LEETEN 01/02/2013-FROM:
			ASSERT_NETCDF(
				nc_get_vara_double(iNcId, iNcVarId, puStarts, puCounts, pdSums) );
			#else	// MOD-BY-LEETEN 01/02/2013-TO:
			ASSERT_NETCDF(
				nc_get_vara(iNcId, ncVarSAT, puStarts, puCounts, (void*)&pdSums[0]) );
			#endif	// MOD-BY-LEETEN 01/02/2013-END

			for(size_t b = 0; b < UGetNrOfBins(); b++)
				vdSums[b] = pdSums[b];

			delete [] pdSums;
		}		
		
		CSATFileEncoder()
		{
		}

		virtual	// ADD-BY-LEETEN 01/02/2013
		~CSATFileEncoder()
		{
			#if	0	// DEL-BY-LEETEN 01/02/2013-BEGIN
			ASSERT_NETCDF(nc_close(iNcId));
			for(size_t d = 0; d < UGetNrOfDims() + 1; d++)
				free(pszNcDimNames[d]);
			#endif	// DEL-BY-LEETEN 01/02/2013-END
		}
	};
}
