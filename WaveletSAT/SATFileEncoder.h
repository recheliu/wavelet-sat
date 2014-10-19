#pragma once

#if	WITH_BOOST
#include <boost/filesystem/operations.hpp>
namespace fs = boost::filesystem;
#endif	// #if	WITH_BOOST

#include <iostream>
#include <vector>
#include <unordered_map>
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
#include "SATFileNetCDF.h"	

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
	template<
		typename DT,				//!< Type of the data
		typename ST = typeSum,		//!< Type of the sum
		typename BT = typeBin		//!< Type of the bin
	>
	class CSATFileEncoder:
		virtual public CHeaderBase,
		virtual public CSATFileNetCDF,	
		virtual public CEncoderBase<DT, ST, BT>
	{
protected:	
		//! The storage to store the original data.
		vector< unordered_map<BT, ST> > vmapHists;

		char szNetCdfPathFilename[NC_MAX_NAME];

		//! The base of all scanlines for all dimensions
		vector< vector<size_t> > vvuSliceScanlineBase;

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
			const BT& uBin, 
			const ST& weight,
			void *_Reserved = NULL
		)
		{
			size_t uIndex = UConvertSubToIndex(vuPos, vuDimLengths);
			unordered_map<BT, ST>& mapHist = vmapHists[uIndex];
			typename unordered_map<BT, ST>::iterator imapHist = mapHist.find(uBin);
			if(mapHist.end() == imapHist )
				mapHist.insert(pair<BT, ST>(uBin, weight));
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
				   (int)(uNrOfDims + 1),
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
			// finish the definition mode
			ASSERT_NETCDF(nc_enddef(iNcId));

			/////////////////////////////////////////
			// start to create the SAT
			vector<size_t> vuScanLineBase;	// the base coordinate of the scane line
			vuScanLineBase.resize(this->UGetNrOfDims());

			vector<size_t> vuOtherDimLengths;	// the base coordinate of the scane line
			vuOtherDimLengths.resize(this->UGetNrOfDims());

			vector<double> pSAT;
			pSAT.assign(this->uDataSize, 0);
			for(size_t 	b = 0;	
					b < UGetNrOfBins(); 
					b++)
			{
				static int iPrintTiming;
				LIBCLOCK_INIT(iPrintTiming, __FUNCTION__);	
				LIBCLOCK_BEGIN(iPrintTiming);	
				for(size_t i = 0; i < this->uDataSize; i++)
				{
					const unordered_map<BT, ST>& mapHist = vmapHists[i];
					typename unordered_map<BT, ST>::const_iterator imapHist = mapHist.find((BT)b);
					pSAT[i] = (double)( mapHist.end() == imapHist )?0:imapHist->second;
				}
				
				LIBCLOCK_END(iPrintTiming);	
				LIBCLOCK_BEGIN(iPrintTiming);	
				for(size_t uOffset = 1, d = 0; d < this->UGetNrOfDims(); uOffset *= vuDimLengths[d], d++)
				{
					size_t uNrOfScanLines = this->uDataSize / this->vuDimLengths[d];
					vuOtherDimLengths = vuDimLengths;
					vuOtherDimLengths[d] = 1;

					for(size_t i = 0; i < uNrOfScanLines; i++)
					{
						size_t uScanLineBase = vvuSliceScanlineBase[d][i];
						for(size_t 	j = 1, uIndex = uScanLineBase; 
								j < vuDimLengths[d]; 
								j++, uIndex += uOffset)
							pSAT[uIndex + uOffset] += pSAT[uIndex];
					}
				}
				LIBCLOCK_END(iPrintTiming);	
				LIBCLOCK_BEGIN(iPrintTiming);	
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

				ASSERT_NETCDF(nc_put_vara(
						   iNcId,
						   ncVarSAT,
						   puStart,
						   puCount,
						   (void*)&pSAT[0]));
				LIBCLOCK_END(iPrintTiming);	
				LIBCLOCK_PRINT(iPrintTiming);	
			}
			ASSERT_NETCDF(nc_close(iNcId));

			iNcId = 0;

		}

		//! Save the coefficients to a file
		virtual 
		void
		_SaveFile
		(
			const char* szFilepathPrefix,
			void *_Reserved = NULL
		)
		{
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
		}

		virtual 
		void
		_DefineBinVar
		(
			const int iNcId,
			const char *szVar,
			const int nctype,
			const int iNrOfDims,
			const int piDimIds[],
			int& ncVar
		)
		{
			ASSERT_NETCDF(nc_def_var(
					iNcId,
					szVar,
					nctype,
					iNrOfDims,
					piDimIds,
					&ncVar));
			#if WITH_NETCDF4
			ASSERT_NETCDF(nc_def_var_deflate(
				   iNcId,
				   ncVar, 
				   0, 
				   1, 
				   iDeflateLevel));
			#endif	// #if WITH_NETCDF4
		}

		virtual 
		void
		_SaveBins
		(
			const char* szFilepathPrefix,
			void *_Reserved = NULL
		)
		{
			const char* szExt;
			#if WITH_NETCDF4
			szExt = "nc4";
			#else	// #if WITH_NETCDF4
			szExt = "nc";
			#endif	// #if WITH_NETCDF4
			char szBinFilepath[NC_MAX_NAME];
			sprintf(szBinFilepath, "%s.bin.%s", szFilepathPrefix, szExt);

			// Create the file for the bins.
			#if !WITH_NETCDF4 
			ASSERT_NETCDF(nc_create(
    				szBinFilepath,
    				NC_CLOBBER,
    				&iNcId));
			#else	// #if !WITH_NETCDF4 
			ASSERT_NETCDF(nc_create(
    				szBinFilepath,
				NC_CLOBBER | NC_NETCDF4,
    				&iNcId));
			#endif	// #if !WITH_NETCDF4 

			// define the dimension for #bins
			ASSERT_NETCDF(nc_def_dim(
						iNcId,
						szDimBin,
						(int)this->UGetNrOfBins(),
						&ncDimBin) );

			// define the dimension of #dimension.
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

			// scan through all maps to decide the total #of bins to save
			size_t uTotalNrOfBins = 0;
			for(size_t v = 0; v < this->vmapHists.size(); v++)
				uTotalNrOfBins += this->vmapHists[v].size();

			// define the dimension for #bins
			ASSERT_NETCDF(nc_def_dim(
						iNcId,
						szDimCoef,
						(int)uTotalNrOfBins,
						&ncDimCoef) );

			// now define the variable for the bins.
			int piNcDimIds[NC_MAX_DIMS];
			for(size_t d = 0; d < uNrOfDims; d++)
				piNcDimIds[d] = vncDimData[uNrOfDims - 1 - d];
			_DefineBinVar
			(
				iNcId,
				this->szVarCoefCount,
				nctypeCoefCount,
				(int)uNrOfDims,
				piNcDimIds,
				ncVarCoefCount
			);
			_DefineBinVar
			(
				iNcId,
				this->szVarCoefOffset,
				nctypeCoefOffset,
				(int)uNrOfDims,
				piNcDimIds,
				ncVarCoefOffset
			);
			_DefineBinVar
			(
				iNcId,
				this->szVarCoefBin,
				nctypeCoefBin,
				1,
				&ncDimCoef,
				ncVarCoefBin
			);
			_DefineBinVar
			(
				iNcId,
				this->szVarCoefSum,
				nctypeCoefSum,
				1,
				&ncDimCoef,
				ncVarCoefSum
			);
			// finish the definition mode
			ASSERT_NETCDF(nc_enddef(iNcId));

			vector<typeCoefCount>	vusCounts;	vusCounts.assign(this->uDataSize,	(typeCoefCount)0);
			vector<typeCoefOffset>	vuOffsets;	vuOffsets.assign(this->uDataSize,	(typeCoefOffset)0);
			vector<typeBin>			vusBins;	vusBins.assign(uTotalNrOfBins,		(typeBin)0);
			vector<typeSum>			vdSums;		vdSums.assign(uTotalNrOfBins,		(typeSum)0);

			for(size_t uOffset = 0, k = 0, i = 0; i < this->uDataSize; uOffset += vmapHists[i].size(), i++)
			{
				vusCounts[i] =	(typeCoefCount)vmapHists[i].size();
				vuOffsets[i] =	(typeCoefOffset)uOffset;
				for(unordered_map<BT, ST>::iterator 
						imap = vmapHists[i].begin();
					imap != vmapHists[i].end();
					imap++, k++)
				{
					vusBins[k]	= (typeBin)imap->first;
					vdSums[k] = (typeSum)imap->second;
				}
			}
			ASSERT_NETCDF(nc_put_var(
						iNcId,
						ncVarCoefCount,
						(void*)&vusCounts.data()[0]));
			ASSERT_NETCDF(nc_put_var(
						iNcId,
						ncVarCoefOffset,
						(void*)&vuOffsets.data()[0]));
			ASSERT_NETCDF(nc_put_var(
						iNcId,
						ncVarCoefBin,
						(void*)&vusBins.data()[0]));
			ASSERT_NETCDF(nc_put_var(
						iNcId,
						ncVarCoefSum,
						(void*)&vdSums.data()[0]));
			ASSERT_NETCDF(nc_close(iNcId));
			iNcId = 0;
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
			vmapHists.resize(uDataSize);

			vvuSliceScanlineBase.resize(UGetNrOfDims());
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				vector<size_t> vuSliceLengths = vuDimLengths;
				vuSliceLengths[d] = 1;
				size_t uNrOfScanlines = uDataSize / vuDimLengths[d];
				vvuSliceScanlineBase[d].resize(uNrOfScanlines);
				for(size_t s = 0; s < uNrOfScanlines; s++)
				{
					vector<size_t> vuScanlineBase;
					_ConvertIndexToSub(s, vuScanlineBase, vuSliceLengths);
					vvuSliceScanlineBase[d][s] = UConvertSubToIndex(vuScanlineBase, vuDimLengths);
				}
			}
		}
		
		//! Compute and display statistics for the computed wavelet coefficients.
		virtual 
		void
		_ShowStatistics
		(
			void *_Reserved = NULL
		)
		{
			#if		WITH_BOOST
			size_t uFileSize = fs::file_size( szNetCdfPathFilename );
			#else	// #if WITH_BOOST
		  FILE *fp;
		  fp = fopen(szNetCdfPathFilename, "rb");
		  fseek(fp, 0, SEEK_END);
		  size_t uFileSize = ftell(fp);
		  fclose(fp);
			#endif	// #if WITH_BOOST	
		  LOG_VAR(uFileSize);
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
				nc_get_vara(iNcId, ncVarSAT, puStarts, puCounts, (void*)&pdSums[0]) );

			for(size_t b = 0; b < UGetNrOfBins(); b++)
				vdSums[b] = pdSums[b];

			delete [] pdSums;
		}		
		
		CSATFileEncoder()
		{
		}

		virtual	
		~CSATFileEncoder()
		{
		}
	};
}

/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
