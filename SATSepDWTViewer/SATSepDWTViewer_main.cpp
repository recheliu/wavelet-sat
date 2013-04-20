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
CSATSepDWT3DView cSATSepDWT3DView;	// ADD-BY-LEETEN 02/05/2013

void 
_GlobalFunc(int iWid, unsigned int uiCbId, va_list vaArgs)
{
	// ADD-BY-LEETEN 02/06/2013-BEGIN
	// 	// pass the new created cluster from cSATSepDWTHistView to cSATSepDWT3DView
		// the passed information: bounding box and color
	// vpairBoxColor

	if( cSATSepDWTHistView.IGetId() == iWid )
	{
		switch(uiCbId)	{
		case CGlutWin::CB_MANUAL:	// not a default GlutWin event
		{
			unsigned int uEvent = va_arg(vaArgs, unsigned int);
			switch(uEvent)	{
			case CSATSepDWTHistView::EVENT_PLOTTING_BOX:
			{
				// read the parameters. eg.
				/* 
				int i = va_arg(vaArgs, int);
				double d = va_arg(vaArgs, double);
				*/
				// read #blocks to add
				int iIsPlottingBlocks = va_arg(vaArgs, int);
				vector< pairBlockColor >& vpairBlockColors = *va_arg(vaArgs, vector< pairBlockColor >*);
				cSATSepDWT3DView._SetBlockColors
					(
						(iIsPlottingBlocks)?true:false,
						vpairBlockColors
					);

			} break;
			}
		} break;
		}
	}
	// ADD-BY-LEETEN 02/06/2013-END
	// ADD-BY-LEETEN 03/20/2013-BEGIN
	if( cSATSepDWT3DView.IGetId() == iWid )
	{
		switch(uiCbId)	{
		case CGlutWin::CB_MANUAL:	// not a default GlutWin event
		{
			unsigned int uEvent = va_arg(vaArgs, unsigned int);
			switch(uEvent)	{
			case CSATSepDWT3DView::EVENT_CURSOR_3D:
			{
				int4 i4Location = va_arg(vaArgs, int4);
				cSATSepDWTHistView._SetQueryLocation(i4Location);
			} break;
			}
		} break;
		}
	}
	// ADD-BY-LEETEN 03/20/2013-END
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
	cSATSepDWTHistView._SetInteger(cSATSepDWTHistView.SIZE_OF_FULL_ARRAYS, (long)iSizeOfFullArrays);
	cSATSepDWTHistView._LoadFile(szNcFilePath);

	/*
	cSATSepDWT3DView._SetData(&cSATSepDWT);
	cSATSepDWT3DView.ICreate("2D/3D Viewer");

		cSATSepDWTHistView._SetData(&cSATSepDWT);
	*/
	cSATSepDWTHistView.ICreate("Integral Histogram Wavelet Viewer");

	// ADD-BY-LEETEN 02/05/2013-BEGIN
	ASSERT_OR_LOG(NULL != szVolFilepath, "Missed szVolFilepath");
	cSATSepDWT3DView._SetData((CSATSepDWT*)&cSATSepDWTHistView);

	cSATSepDWT3DView._LoadData(szVolFilepath);
	cSATSepDWT3DView.ICreate("3D Viewer");
	// ADD-BY-LEETEN 02/05/2013-END

	CGlutWin::_RegisterGlobalFunc(_GlobalFunc);	// ADD-BY-LEETEN 02/06/2013

	glutMainLoop();

	return 0;
}

