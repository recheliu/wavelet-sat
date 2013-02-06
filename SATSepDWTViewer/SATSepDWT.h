#include <string>
using namespace std;

#include <GL/glew.h>	

#include <vector_types.h>	// ADD-BY-LEETEN 02/06/2013

#include "libopt.h"
#include "liblog.h"
#include "libclock.h"	

#include "SimpleNDFile.h"

typedef CSimpleNDFile<double> CSATSepDWT;

// ADD-BY-LEETEN 02/06/2013-BEGIN
//! The type that pairs a box (of two int4 as the corners) and colors
typedef pair< pair<float4, float4>, float4 > pairBlockColor;
// ADD-BY-LEETEN 02/06/2013-END