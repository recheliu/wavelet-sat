#include <GL/glew.h>

#include "shader.h"	
#include "libbuf3d.h"

#include "libclock.h"	

#include "SATSepDWTQueryView.h"
using namespace WaveletSAT;
using namespace SATSepDWT;

//////////////////// CGlutWin methods //////////////////// 
void 
CQueryView::_IdleFunc()
{
}

void 
CQueryView::
	_MotionFunc
	(
		int x, 
		int y
	)
{
}

void 
CQueryView::_MouseFunc
	(
		int button, 
		int state, 
		int x, 
		int y
	)
{
}

void 
CQueryView::
	_InitGl()
{
}

void 
CQueryView::_InitFunc()
{
	_KeepUpdateOn();
	
	/////////////////////////////////////////////
	// set up GLUI
	GLUI *pcGlui = PCGetGluiSubwin();

	cQuery._AddGlui(this, pcGlui, NULL, pcBlockTree->VGetDimLengths());

}

void 
CQueryView::
	_TimerFunc(unsigned short usEvent)
{
}

void 
CQueryView::_DisplayFunc()
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
	_DrawString(szHistogramMax);		

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
CQueryView::_KeyboardFunc(unsigned char ubKey, int x, int y)
{
}

void 
CQueryView::_ReshapeFunc(int w, int h)
{
}

void 
CQueryView::_GluiFunc(unsigned short usValue)
{
	switch(usValue)
	{
	case GLUI_EVENT_QUERY:
	{
		_QueryHistogram();
	} break;
	}
}

void 
CQueryView::
_QueryHistogram
(
	void *_Reserved
)
{
	LIBCLOCK_INIT(cQuery.iIsPrintingTiming, "GLUI_EVENT_QUERY");
	LIBCLOCK_BEGIN(cQuery.iIsPrintingTiming);
	vector<size_t> vCorner0, vCorner1;
	vCorner0.resize(pcBlockTree->UGetNrOfDims());
	vCorner1.resize(pcBlockTree->UGetNrOfDims());
	int *piLocation =	&cQuery.pcQuery->i4Location.x;
	int *piSize =		&cQuery.pcQuery->i4Size.x;
	for(size_t d = 0; d < pcBlockTree->UGetNrOfDims(); d++){
		vCorner0[d] = min(max(piLocation[d], 0), (int)pcBlockTree->VGetDimLengths()[d] - 1);
		vCorner1[d] = min(max(piLocation[d] + piSize[d], 0), (int)pcBlockTree->VGetDimLengths()[d] - 1);
	}

	vector<WaveletSAT::typeSum> vdH;
	vdH.resize(pcBlockTree->UGetNrOfBins());
	pcBlockTree->_GetRegionSums(vCorner0, vCorner1, vdH);
	LIBCLOCK_END(cQuery.iIsPrintingTiming);

	LIBCLOCK_BEGIN(cQuery.iIsPrintingTiming);
	_SetHistgroam(vdH);
	LIBCLOCK_END(cQuery.iIsPrintingTiming);
	LIBCLOCK_PRINT(cQuery.iIsPrintingTiming);
}

//////////////////// Constructors/Destructors //////////////////// 
CQueryView::
	CQueryView(void)
{
	// add a panel for the UI control
	_AddGluiSubwin(GLUI_SUBWINDOW_TOP);
}

CQueryView::
	~CQueryView(void)
{
}

