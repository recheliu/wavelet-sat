#pragma once

#include <string>
using namespace std;

#include <GL/glew.h>	
#include <GL/glut.h>

#include <vector_types.h>	
#include <vector_functions.h>	

#include "libopt.h"
#include "liblog.h"
#include "libclock.h"	

#include "SimpleNDFile.h"

namespace SATSepDWT 
{
	struct CQuery
	{
		int iIsPlot;
		int4 i4Location;
		int4 i4Size;
	};
}