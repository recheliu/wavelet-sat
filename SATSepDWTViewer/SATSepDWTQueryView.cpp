#include <GL/glew.h>

#include "shader.h"	
#include "libbuf3d.h"

#include "libclock.h"	

#include "SATSepDWTQueryView.h"

using namespace WaveletSAT;

//////////////////// CGlutWin methods //////////////////// 
void 
CSATSepDWTQueryView::_IdleFunc()
{
}

void 
CSATSepDWTQueryView::
	_MotionFunc
	(
		int x, 
		int y
	)
{
}

void 
CSATSepDWTQueryView::_MouseFunc
	(
		int button, 
		int state, 
		int x, 
		int y
	)
{
}

void 
CSATSepDWTQueryView::
	_InitGl()
{
}

void 
CSATSepDWTQueryView::_InitFunc()
{
	_KeepUpdateOn();
	
	/////////////////////////////////////////////
	// set up GLUI
	GLUI *pcGlui = PCGetGluiSubwin();

}

void 
CSATSepDWTQueryView::
	_TimerFunc(unsigned short usEvent)
{
}

void 
CSATSepDWTQueryView::_DisplayFunc()
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if( vdHistogram.empty() )
		return;

	// convert the screen space to [0 ... #Bins] x [0 ... #Levels]
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	char szHistogramMax[1024];
	sprintf(szHistogramMax, "%.2f", dHistogramMax);
	_DrawString3D(szHistogramMax, -1.0f, 0.9f);		

	glTranslatef(-1.0f, -1.0f, 0.0f);
	glScalef(2.0f/(float)vdHistogram.size(), 2.0f/(float)dHistogramMax, 1.0f);

	glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
	glBegin(GL_LINE_STRIP);
	for(size_t b = 0; b < vdHistogram.size(); b++)
	{
		glVertex2f((float)b, vdHistogram[b]);
		glVertex2f((float)(b + 1), vdHistogram[b]);
	}
	glEnd();

	// restore the coordinate system back
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void 
CSATSepDWTQueryView::_KeyboardFunc(unsigned char ubKey, int x, int y)
{
}

void 
CSATSepDWTQueryView::_ReshapeFunc(int w, int h)
{
}

void 
CSATSepDWTQueryView::_GluiFunc(unsigned short usValue)
{
	switch(usValue)
	{
	}
}

//////////////////// Constructors/Destructors //////////////////// 
CSATSepDWTQueryView::
	CSATSepDWTQueryView(void)
{
	// add a panel for the UI control
	_AddGluiSubwin(GLUI_SUBWINDOW_LEFT);
}

CSATSepDWTQueryView::
	~CSATSepDWTQueryView(void)
{
}

