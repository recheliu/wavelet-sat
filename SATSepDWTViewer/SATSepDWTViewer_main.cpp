#include <string>
using namespace std;

#include <GL/glew.h>	

#include "libopt.h"
#include "liblog.h"
#include "libclock.h"	

#include "SATSepDWT.h"
#include "SATSepDWT3DView.h"
#include "SATSepDWTHistView.h"

CSATSepDWTHistView cSATSepDWTHistView;
// CSATSepDWT3DView cSATSepDWT3DView;

void
_UpdateTf()
{
}

void 
_GlobalFunc(int iWid, unsigned int uiCbId, va_list vaArgs)
{
}

int
main(int argn, char* argv[])
{
	// Initialize GlutWin first
	CGlutWin::_Init(
		&argn, argv, 
		GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_STENCIL);

	// parse the arguments
	_OPTInit(true);

	char* szNcFilePath = NULL;
	_OPTAddStringVector(
		"--nc-filepath", 1,
		&szNcFilePath, szNcFilePath);

	int iSizeOfFullArrays = 0;
	_OPTAddIntegerVector(
		"--size-of-full-arrays", 1,
		&iSizeOfFullArrays, iSizeOfFullArrays);
	_OPTAddComment("--size-of-full-arrays", 
		"Size (in MB) of the full arrays from all bin SATs");

	bool bIsOptParsed = BOPTParse(argv, argn, 1);
	ASSERT_OR_LOG(bIsOptParsed, "Invalid options.");
	ASSERT_OR_LOG(szNcFilePath, "THe input file is missed.");

	// load the data
	cSATSepDWTHistView._SetInteger(cSATSepDWTHistView.SIZE_OF_FULL_ARRAYS, (long)iSizeOfFullArrays);
	cSATSepDWTHistView._LoadFile(szNcFilePath);

	/*
	cSATSepDWT3DView._SetData(&cSATSepDWT);
	cSATSepDWT3DView.ICreate("2D/3D Viewer");

		cSATSepDWTHistView._SetData(&cSATSepDWT);
	*/
	cSATSepDWTHistView.ICreate("Integral Histogram Wavelet Viewer");

	glutMainLoop();

	return 0;
}

