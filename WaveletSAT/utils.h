#pragma once

// ADD-BY-LEETEN 01/31/2013-BEGIN
#include <type_traits>
using namespace std;
// ADD-BY-LEETEN 01/31/2013-END

#define _USE_MATH_DEFINES
#include <math.h>

#if defined (WIN32)
	#include <psapi.h>	
	#pragma comment (lib, "psapi.lib")
#else	// #if defined (WIN32)
	#include <sys/time.h>
	#include <sys/resource.h>
#endif	// #if defined (WIN32)

// ADD-BY-LEETEN 01/23/2013-BEGIN
#include <netcdf.h>
#include "NrrdIO.h"

#include "gnuplot.h"
// ADD-BY-LEETEN 01/23/2013-END

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	// ADD-BY-LEETEN 01/23/2013-BEGIN
	namespace Statistics 
	{
		template<typename ST,	typename BT, typename WT>
		ST 
		Mean(
			const vector< pair< BT, WT > >& vpairBins
		)
		{
			ST Sum = (ST)0;
			ST Mean = (ST)0;
			for(typename vector< pair< BT, WT > >::const_iterator 
					ivpairBin = vpairBins.begin();
				ivpairBin != vpairBins.end();
				ivpairBin++)
			{
				ST Bin = (ST)ivpairBin->first;
				ST Value = (ST)ivpairBin->second;
				Mean += Bin * Value;
				Sum += Value;
			}
			return Mean/Sum;
		}

		template<typename ST, typename BT, typename WT>
		ST 
		Entropy(
			const vector< pair< BT, WT > >& vpairBins
		)
		{
			ST S = (ST)0;
			ST E = (ST)0;
			for(typename vector< pair< BT, WT > >::const_iterator 
					ivpairBin = vpairBins.begin();
				ivpairBin != vpairBins.end();
				ivpairBin++)
			{
				ST Value = (ST)ivpairBin->second;
				S += Value;
				E += Value * (ST)log((double)Value);
			}
			E = -E / S + (ST)log((double)S);
			E /= M_LN2;

			return E;
		}

		template<typename ST,	typename BT, typename WT>
		ST 
		Count(
			const vector< pair< BT, WT > >& vpairBins
		)
		{
			return (ST)vpairBins.size();
		}

		template<typename ST,	typename BT, typename WT>
		ST 
		StdDev(
			const vector< pair< BT, WT > >& vpairBins
		)
		{
			ST SquaredSum = (ST)0;
			ST Mean = (ST)0;
			ST Sum = (ST)0;
			for(typename vector< pair< BT, WT > >::const_iterator 
					ivpairBin = vpairBins.begin();
				ivpairBin != vpairBins.end();
				ivpairBin++)
			{
				ST Bin = (ST)ivpairBin->first;
				ST Value = (ST)ivpairBin->second;
				Mean += Bin * Value;
				SquaredSum += Bin * Bin * Value;
				Sum += Value;
			}
			SquaredSum /= Sum;
			Mean /= Sum;
			ST Var = SquaredSum - Mean * Mean;
			return (ST)sqrt((double)Var);
		}
	};

	template<typename NT>
	void
	_SaveNrrd
	(
		const vector<size_t>& vuDimLengths,
		void* pData,
		const char* szFilepath,
		void *_Reserved = NULL
		)
	{
		int nrrdType;
		if( is_same<NT, double>::value )
			nrrdType = nrrdTypeDouble;
		if( is_same<NT, float>::value )
			nrrdType = nrrdTypeFloat;

		// create the .nhdr filename
		char pathNhdr[MAX_NC_NAME+1];      
		strcpy(pathNhdr, szFilepath);
		char *szExtension = strrchr(pathNhdr, '.');
		if( 0 != strcmp(szExtension, ".nhdr" ) )
			strcat(pathNhdr, ".nhdr");

		// generate the corresponding .raw filename
		char pathRaw[MAX_NC_NAME+1];       
		strcpy(pathRaw, pathNhdr);      
		strcpy(strrchr(pathRaw, '.'), ".raw");

		// generate the base filename for .nhdr
		char *szFilename = strrchr(pathNhdr, '/');
		#if defined(WIN32)
		if( !szFilename )
		szFilename = strrchr(pathNhdr, '\\');
		#endif // #if defined(WIN32)

		szFilename = ( szFilename )?(szFilename+1):&pathNhdr[0];
		char pathLeafNhdr[MAX_NC_NAME];  
		strcpy(pathLeafNhdr, szFilename);
		const char *szLocalNhdr = &pathLeafNhdr[0];

		// generate the base filename for .raw
		char pathLeafRaw[MAX_NC_NAME];   
		strcpy(pathLeafRaw, pathLeafNhdr);
		strcpy(strrchr(pathLeafRaw, '.'), ".raw");

		size_t uDataLength = 1;
		for(size_t d = 0; d < vuDimLengths.size(); d++)
			uDataLength *= vuDimLengths[d];

		Nrrd *nrrdOut = nrrdNew();
		nrrdWrap_nva(nrrdOut, pData, nrrdType, (unsigned int)vuDimLengths.size(), vuDimLengths.data());
		nrrdSave(szLocalNhdr, nrrdOut, NULL);

		// nrrdIoStateNix(nioOut);
		nrrdNix(nrrdOut);

		if( 0 != strcmp( pathLeafNhdr, pathNhdr ) )
		{
			remove(pathNhdr);
			remove(pathRaw);
			rename(pathLeafNhdr, pathNhdr);
			rename(pathLeafRaw, pathRaw);
		} 

		#if	WITH_GNUPLOT
		// now apply GNUPLOT
		char szPngPathFilename[MAX_NC_NAME+1];      
		strcpy(szPngPathFilename, pathNhdr);
		strcat(szPngPathFilename, ".png");

		// launch GNUPLOT to plot the FTLE
		static char szGnuplotCommand[2048];
		sprintf(szGnuplotCommand, "%s -e \""
			// ADD-BY-LEETEN 10/06/2011-BEGIN
			"set autoscale;"
			"unset label;"
			"set terminal png;"		
			"set output '%s';"		
			"set size ratio %f;"	
			"set xrange [0:%d];"
			"set yrange [0:%d];"
			"plot '%s' binary array=(%d,%d) format='%%lf' with image;"
			, 
			GNUPLOT_EXE_PATHFILENAME, 
			szPngPathFilename,	
			(float)vuDimLengths[1]/(float)vuDimLengths[0],	
			vuDimLengths[0] - 1, 
			vuDimLengths[1] - 1,	
			pathRaw, vuDimLengths[0],	vuDimLengths[1]
			);
		system(szGnuplotCommand);
		#endif	// #if	WITH_GNUPLOT
	}
	// ADD-BY-LEETEN 01/23/2013-END

	// ADD-BY-LEETEN 10/30/2012-BEGIN
	inline	// ADD-BY-LEETEN 01/27/2012
	void
	_ConvertIndexToSub
	(
		size_t uIndex,
		vector<size_t>& vuSub,
		const vector<size_t>& vuDimLengths,
		void* _Reserved = NULL
	)
	{
		size_t uNrOfDims = vuDimLengths.size();
		if( uNrOfDims != vuSub.size() )
			vuSub.resize(uNrOfDims);
		for(size_t d = 0; d < uNrOfDims; d++)
		{
			vuSub[d] = uIndex % vuDimLengths[d];
			uIndex /= vuDimLengths[d];
		}
	};

	// ADD-BY-LEETEN 2013/12/07-BEGIN
	inline 
	size_t
	UGetProduct
	(
		const vector<size_t> vuLengths,
		void* _Reserved = NULL
	)
	{
		size_t uProd = 1;
		for(size_t d = 0; d < vuLengths.size(); d++)
			uProd *= vuLengths[d];
		return uProd;
	}
	// ADD-BY-LEETEN 2013/12/07-END
	inline	// ADD-BY-LEETEN 01/27/2012
	size_t
	UConvertSubToIndex
	(
		const vector<size_t>& vuSub,
		const vector<size_t>& vuDimLengths,
		void* _Reserved = NULL
	)
	{
		size_t uNrOfDims = vuDimLengths.size();
		size_t uIndex = 0;
		for(size_t d = 0, vuSubSize = 1; d < vuSub.size(); vuSubSize *= vuDimLengths[d], d++)
			uIndex += vuSubSize * vuSub[d];
		return uIndex;
	}
	// ADD-BY-LEETEN 10/30/2012-END

inline
size_t
UGetMemoryUsage() 
{
	#if defined(WIN32)
	  PROCESS_MEMORY_COUNTERS memCounter;
	  BOOL result = GetProcessMemoryInfo(
					     GetCurrentProcess(),
					     &memCounter,	
					     sizeof( memCounter ));
	  return memCounter.WorkingSetSize;
	#else		// #if defined(WIN32)
	  int who = RUSAGE_SELF;
	  struct rusage usage;	
	  int ret;			   
	  getrusage(who,&usage);
	  return usage.ru_maxrss;
	#endif		// #if defined(WIN32)
}

#define _ShowMemoryUsage(bIsOutputToError) \
	{				   \
		size_t uMemoryUsage = UGetMemoryUsage();	\
	  if( bIsOutputToError )			\
	    {						\
	      LOG_VAR_TO_ERROR(uMemoryUsage);	\
	    }						\
	  else						\
	    LOG_VAR(uMemoryUsage);			\
	}

}

/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
