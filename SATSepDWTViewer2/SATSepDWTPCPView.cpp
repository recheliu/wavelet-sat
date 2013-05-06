#include <algorithm>
#include <GL/glew.h>

#include "shader.h"	
#include "libbuf3d.h"

#include "libclock.h"	

#include "SATSepDWTPCPView.h"
using namespace SATSepDWT;

//////////////////// CGlutWin methods //////////////////// 
void 
CPCPView::_InitFunc()
{
	_KeepUpdateOn();
	
	/////////////////////////////////////////////
	// set up GLUI
	GLUI *pcGlui = PCGetGluiSubwin();

	GLUI_Panel *pcPanel_Line = pcGlui->add_rollout("Line");
	GLUI_Panel *pcPanel_Color = pcGlui->add_panel_to_panel(pcPanel_Line, "Color");
	static char* pszChannels[] = {"R", "G", "B", "A"};
	float *pfColor = &f4Color.x;
	for(int c = 0; c < sizeof(pszChannels)/sizeof(pszChannels[0]); c++)
	{
		GLUI_Spinner* pcSpinner = pcGlui->add_spinner_to_panel(pcPanel_Color, pszChannels[c], GLUI_SPINNER_FLOAT, &pfColor[c]);
		pcSpinner->set_float_limits(0.0f, 1.0f);
	}
	GLUI_Spinner* pcSpinner_Width = pcGlui->add_spinner_to_panel(pcPanel_Line, "Width", GLUI_SPINNER_FLOAT, &fWidth);
		pcSpinner_Width->set_float_limits(1.0f, 16.0f);

	cFilter._AddGlui(this, pcGlui, NULL, pcBlockTree->uMaxLevel);
}

void 
CPCPView::_DisplayFunc()
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// convert the screen space to [0 ... #Bins] x [0 ... #Levels]
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	// scale the coordinate sysm s.t. the range of X axis is [-1, uMaxLevel] 
	size_t uNrOfLevels = pcBlockTree->uMaxLevel;
	glTranslatef(-1.0f, -1.0f, 0.0f);
	glScalef(2.0f/(float)(uNrOfLevels + 1), 2.0f, 1.0f);
	glTranslatef(+1.0f, 0.0f, 0.0f);

	// plot the axis
	glBegin(GL_LINES);
	glColor4fv(&f4Color.x);
	for(size_t l = 0; l < uNrOfLevels; l++)
	{
		glVertex2f((float)l, 0.0f);
		glVertex2f((float)l, 1.0f);
	}
	glEnd();

	glPushAttrib(
		GL_COLOR_BUFFER_BIT |
		0 );
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);

	pcBlockTree->_RenderPolylines(CBlock::MODE_NONE, this->f4Color);
	pcBlockTree->_RenderPolylines(CBlock::MODE_ASSIGNED, this->f4Color);

	glPopAttrib();	
	//	GL_COLOR_BUFFER_BIT |

	// restore the coordinate system back
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void 
CPCPView::_GluiFunc(unsigned short usValue)
{
	switch(usValue)
	{
	case GLUI_EVENT_FILTER_ASSIGN:
	{
	} break;
	}
}

//////////////////// Constructors/Destructors //////////////////// 
CPCPView::
	CPCPView(void)
{
	f4Color = make_float4(0.0f, 0.0f, 0.0f, 1.0f);

	// add a panel for the UI control
	_AddGluiSubwin(GLUI_SUBWINDOW_LEFT);
}

CPCPView::
	~CPCPView(void)
{
}

