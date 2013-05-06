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

#include "SATSepDWTBlock.h";
using namespace SATSepDWT;
CSimpleNDFile<CBlock::typeData> *CBlock::pcSimpleND;
