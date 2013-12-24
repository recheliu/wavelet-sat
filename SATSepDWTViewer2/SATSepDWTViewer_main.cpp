#include <string>
using namespace std;

#include <GL/glew.h>	

#include "libopt.h"
#include "liblog.h"
#include "libclock.h"	

// #include "SATSepDWT.h"
#include "SATSepDWTBlockTree.h"
#include "SATSepDWTQuery.h"
#include "SATSepDWT3DView.h"
#include "SATSepDWTHistView.h"
#include "SATSepDWTPCPView.h"
#include "SATSepDWTQueryView.h"
using namespace SATSepDWT;

CQuery cQuery;
CBlockTree cBlockTree;
CHistView cHistView;
CPCPView cPCPView;
CQueryView cQueryView;
C3DView c3DView;

void 
_GlobalFunc(int iWid, unsigned int uiCbId, va_list vaArgs)
{
	if( c3DView.IGetId() == iWid )
	{
		switch(uiCbId)	{
		case CGlutWin::CB_MANUAL:	// not a default GlutWin event
		{
			unsigned int uEvent = va_arg(vaArgs, unsigned int);
			switch(uEvent)	{
			case C3DView::EVENT_CURSOR_3D:
			{
				cQueryView._QueryHistogram();
			} break;
			}
		} break;
		}
	}
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

	// ADD-BY-LEETEN 02/05/2013-BEGIN
	char* szVolFilepath = NULL;
	_OPTAddStringVector(
		"--vol-filepath", 1,
		&szVolFilepath, szVolFilepath);
	// ADD-BY-LEETEN 02/05/2013-END

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
	cBlockTree._SetInteger(cBlockTree.SIZE_OF_FULL_ARRAYS, (long)iSizeOfFullArrays);
	cBlockTree._LoadFile(szNcFilePath);

	cHistView._SetBlockTree(&cBlockTree);
	cHistView.ICreate("Integral Histogram Wavelet Viewer");

	cPCPView._SetBlockTree(&cBlockTree);
	cPCPView.ICreate("PCP Viewer");

	cQueryView._SetBlockTree(&cBlockTree);
	cQueryView._SetQuery(&cQuery);
	cQueryView.ICreate("Query Viewer");

	c3DView._SetBlockTree(&cBlockTree);
	c3DView._SetQuery(&cQuery);
	c3DView._LoadData(szVolFilepath);
	c3DView.ICreate("3D Viewer");

	CGlutWin::_RegisterGlobalFunc(_GlobalFunc);	

	glutMainLoop();

	return 0;
}


/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
