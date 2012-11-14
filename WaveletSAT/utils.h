#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

#if defined (WIN32)
	#include <psapi.h>	
	#pragma comment (lib, "psapi.lib")
#else	// #if defined (WIN32)
	#include <sys/time.h>
	#include <sys/resource.h>
#endif	// #if defined (WIN32)

/*
Usage: The application just calls _SetDimLengths() first and then _AllocateBins() to setup the class. 
Then the user call _Update(vuPos, value) to update the value at the given location.
*/
namespace WaveletSAT
{
	// ADD-BY-LEETEN 10/30/2012-BEGIN
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

#if 0 // MOD-BY-LEETEN 11/13/2012-FROM:
	void
	_ShowMemoryUsage
	(
		bool bIsOutputToError = false,	// ADD-BY-LEETEN 10/30/2012
		void* _Reserved = NULL
	)
	{
		#if defined(WIN32)
		PROCESS_MEMORY_COUNTERS memCounter;
		BOOL result = GetProcessMemoryInfo(
				GetCurrentProcess(),
				&memCounter,
				sizeof( memCounter ));
		// ADD-BY-LEETEN 10/30/2012-BEGIN
		if( bIsOutputToError )
		{
			LOG_VAR_TO_ERROR(memCounter.WorkingSetSize);
		}
		else
		// ADD-BY-LEETEN 10/30/2012-END
		LOG_VAR(memCounter.WorkingSetSize);

		#else	// #if defined(WIN32)
		int who = RUSAGE_SELF; 
		struct rusage usage; 
		int ret; 
		getrusage(who,&usage);
		// ADD-BY-LEETEN 10/30/2012-BEGIN
		if( bIsOutputToError )
		  { // ADD-BY-LEETEN 10/31/2012
			LOG_VAR_TO_ERROR(usage.ru_maxrss);
		  } // ADD-BY-LEETEN 10/31/2012
		else
		// ADD-BY-LEETEN 10/30/2012-END
		LOG_VAR(usage.ru_maxrss);

		#endif	// #if defined(WIN32)
	}
#else // MOD-BY-LEETEN 11/13/2012-TO:

#if defined(WIN32)
#define _ShowMemoryUsage(bIsOutputToError) \
	{								\
	  PROCESS_MEMORY_COUNTERS memCounter;				\
	  BOOL result = GetProcessMemoryInfo(				\
					     GetCurrentProcess(),	\
					     &memCounter,		\
					     sizeof( memCounter ));	\
	  if( bIsOutputToError )					\
	    {								\
	      LOG_VAR_TO_ERROR(memCounter.WorkingSetSize);		\
	    }								\
	  else								\
	    LOG_VAR(memCounter.WorkingSetSize);				\
	}

#else	// #if defined(WIN32)

#define _ShowMemoryUsage(bIsOutputToError) \
	{				   \
	  int who = RUSAGE_SELF;	   \
	  struct rusage usage;		   \
	  int ret;			   \
	  getrusage(who,&usage);	   \
	  if( bIsOutputToError )			\
	    {						\
	      LOG_VAR_TO_ERROR(usage.ru_maxrss);	\
	    }						\
	  else						\
	    LOG_VAR(usage.ru_maxrss);			\
	}

#endif	// #if defined(WIN32)

#endif // MOD-BY-LEETEN 11/13/2012-END
}
