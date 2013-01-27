#include <GL/glew.h>

#include "shader.h"	
#include "libbuf3d.h"

#include "libclock.h"	

#include "SATSepDWT3DView.h"


//////////////////// CGlutWin methods //////////////////// 
void 
CSATSepDWT3DView::_IdleFunc()
{
}

void 
CSATSepDWT3DView::_MouseFunc
	(
		int button, 
		int state, 
		int x, 
		int y
	)
{
}

void 
CSATSepDWT3DView::
	_InitGl()
{
}

void 
CSATSepDWT3DView::_InitFunc()
{
	CDvrWin2::_InitFunc();

	_KeepUpdateOn();
}

void 
CSATSepDWT3DView::
	_TimerFunc(unsigned short usEvent)
{
}

void 
CSATSepDWT3DView::_DisplayFunc()
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);
}

void 
CSATSepDWT3DView::_KeyboardFunc(unsigned char ubKey, int x, int y)
{
}

void 
CSATSepDWT3DView::_ReshapeFunc(int w, int h)
{
	CDvrWin2::_ReshapeFunc(w, h);
	/*
	CClipVolume::_ReshapeFunc(w, h);
	_Redisplay();
	*/
}

void 
CSATSepDWT3DView::_GluiFunc(unsigned short usValue)
{
}

//////////////////////////////////////////////////////
void 
CSATSepDWT3DView::
	_BeginDisplay()
{
}

void 
CSATSepDWT3DView::
	_EndDisplay()
{
}

void 
CSATSepDWT3DView::
	_RenderSlab(
		int iSlab, int iNrOfSlabs,
		double pdModelviewMatrix[], double pdProjectionMatrix[], int piViewport[],
		double dMinX, double dMaxX, 
		double dMinY, double dMaxY, 
		double dMinZ, double dMaxZ
	)
{
}

//////////////////// Constructors/Destructors //////////////////// 
CSATSepDWT3DView::
	CSATSepDWT3DView(void)
{
}

CSATSepDWT3DView::
	~CSATSepDWT3DView(void)
{
}

