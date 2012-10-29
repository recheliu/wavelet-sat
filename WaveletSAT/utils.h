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
	void
	_ShowMemoryUsage
	(
		void* _Reserved = NULL
	)
	{
		#if defined(WIN32)
		PROCESS_MEMORY_COUNTERS memCounter;
		BOOL result = GetProcessMemoryInfo(
				GetCurrentProcess(),
				&memCounter,
				sizeof( memCounter ));
		LOG_VAR(memCounter.WorkingSetSize);

		#else	// #if defined(WIN32)
		int who = RUSAGE_SELF; 
		struct rusage usage; 
		int ret; 
		getrusage(who,&usage);
		LOG_VAR(usage.ru_maxrss);

		#endif	// #if defined(WIN32)
	}
}
