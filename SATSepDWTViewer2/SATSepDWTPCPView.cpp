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

	cFilter._AddGlui(this, pcGlui, NULL, pcBlockTree->uMaxLevel + 1);
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
	size_t uNrOfLevels = pcBlockTree->uMaxLevel + 1;
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

	// ADD-BY-LEETEN 05/07/2013-BEGIN
	// plot the filter boxes
	if( cFilter.iIsActive )
	{
		// ADD-BY-LEETEN 06/23/2013-BEGIN
		CBlock::_ClearBlocksRenderedByPCP();
		// ADD-BY-LEETEN 06/23/2013-END
		glColor4fv(&cFilter.f4Color.x);
		pcBlockTree->_FilteredPaths(
			CBlock::FILTER_ACTION_RENDER,
			cFilter.vf2LevelProb,
			(size_t)cFilter.iTarget,
			cFilter.f4Color);
	}
	// ADD-BY-LEETEN 05/07/2013-END

	glPopAttrib();	
	//	GL_COLOR_BUFFER_BIT |

	// ADD-BY-LEETEN 05/07/2013-BEGIN
	// plot the filter boxes
	if( cFilter.iIsActive )
	{
		glPushAttrib(
			GL_LINE_BIT |
			0);
		glLineStipple(1, 0xcccc);
		glEnable(GL_LINE_STIPPLE);
		glColor4fv(&cFilter.f4Color.x);
		for(size_t l = 0; l < cFilter.vf2LevelProb.size(); l++)
		{
			glBegin(GL_LINE_LOOP);
			glVertex2f((float)l - cFilter.fWidth, cFilter.vf2LevelProb[l].x);
			glVertex2f((float)l + cFilter.fWidth, cFilter.vf2LevelProb[l].x);
			glVertex2f((float)l + cFilter.fWidth, cFilter.vf2LevelProb[l].y);
			glVertex2f((float)l - cFilter.fWidth, cFilter.vf2LevelProb[l].y);
			glEnd();
		}
		glPopAttrib();
			// GL_LINE_BIT |
		glBegin(GL_LINE_LOOP);
		glVertex2f((float)cFilter.iTarget - cFilter.fWidth, cFilter.vf2LevelProb[cFilter.iTarget].x);
		glVertex2f((float)cFilter.iTarget + cFilter.fWidth, cFilter.vf2LevelProb[cFilter.iTarget].x);
		glVertex2f((float)cFilter.iTarget + cFilter.fWidth, cFilter.vf2LevelProb[cFilter.iTarget].y);
		glVertex2f((float)cFilter.iTarget - cFilter.fWidth, cFilter.vf2LevelProb[cFilter.iTarget].y);
		glEnd();
	}
	// ADD-BY-LEETEN 05/07/2013-END
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
		// ADD-BY-LEETEN 05/07/2013-BEGIN
		if( cFilter.iIsActive )
		{
			pcBlockTree->_FilteredPaths(
				CBlock::FILTER_ACTION_ASSIGN,
				cFilter.vf2LevelProb,
				(size_t)cFilter.iTarget,
				cFilter.f4Color
				);
		}
		// ADD-BY-LEETEN 05/07/2013-END
	} break;
	// ADD-BY-LEETEN 05/07/2013-BEGIN
	case GLUI_EVENT_FILTER_RESET:
	{
		if( cFilter.iIsActive )
		{
			pcBlockTree->_FilteredPaths(
				CBlock::FILTER_ACTION_RESET,
				cFilter.vf2LevelProb,
				(size_t)cFilter.iTarget,
				cFilter.f4Color
				);
		}
	} break;
	// ADD-BY-LEETEN 05/07/2013-END
	
	// ADD-BY-LEETEN 06/23/2013-BEGIN
	case GLUI_EVENT_FILTER_ACTIVE:
		if( !this->cFilter.iIsActive )
			CBlock::_ClearBlocksRenderedByPCP();
		break;
	case GLUI_EVENT_FILTER_COLOR:
		CBlock::_SetPCPColor(cFilter.f4Color);
		break;
	// ADD-BY-LEETEN 06/23/2013-END
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


/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
