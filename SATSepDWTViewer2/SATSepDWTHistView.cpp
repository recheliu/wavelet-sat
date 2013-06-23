#include <GL/glew.h>

#include "shader.h"	
#include "libbuf3d.h"

#include "libclock.h"	

#include "SATSepDWTHistView.h"
using namespace SATSepDWT;
using namespace WaveletSAT;	// ADD-BY-LEETEN 03/17/2013

//////////////////// CGlutWin methods //////////////////// 
void 
CHistView::_IdleFunc()
{
}

void 
CHistView::
	_MotionFunc
	(
		int x, 
		int y
	)
{
}

void 
CHistView::_MouseFunc
	(
		int button, 
		int state, 
		int x, 
		int y
	)
{
}

void 
CHistView::
	_InitGl()
{
}

void 
CHistView::_InitFunc()
{
	_KeepUpdateOn();
	
	i2BinRange.x = 0;
	i2BinRange.y = pcBlockTree->UGetNrOfBins() - 1;

	/////////////////////////////////////////////
	// set up GLUI
	GLUI *pcGlui = PCGetGluiSubwin();

	GLUI_Panel *pcPanel_Basic = pcGlui->add_rollout("Basic");
	GLUI_Spinner* pcSpinner_MaxLevel = pcGlui->add_spinner_to_panel(pcPanel_Basic, "Max Level", GLUI_SPINNER_INT, &iNrOfLevelsToDisplay);
		pcSpinner_MaxLevel->set_int_limits(0, (int)pcBlockTree->uMaxLevel);
	GLUI_Spinner* pcSpinner_MaxBin = pcGlui->add_spinner_to_panel(pcPanel_Basic, 
			"Max Bin", GLUI_SPINNER_INT, &i2BinRange.y, 
			IAddWid(GLUI_EVENT_BIN_RANGE), CGlutWin::_GluiCB_static);
		pcSpinner_MaxBin->set_int_limits(0, pcBlockTree->UGetNrOfBins() - 1);
	GLUI_Spinner* pcSpinner_MinBin = pcGlui->add_spinner_to_panel(pcPanel_Basic, 
			"Min Bin", GLUI_SPINNER_INT, &i2BinRange.x,
			IAddWid(GLUI_EVENT_BIN_RANGE), CGlutWin::_GluiCB_static);
		pcSpinner_MinBin->set_int_limits(0,	 pcBlockTree->UGetNrOfBins() - 1);

	{
		GLUI_Panel *pcPanel_Color = pcGlui->add_panel_to_panel(pcPanel_Basic, "Default Color");
		static char* pszChannels[] = {"R", "G", "B", "A"};
		float *pfColor = &f4DefaultColor.x;
		for(int c = 0; c < sizeof(pszChannels)/sizeof(pszChannels[0]); c++)
		{
			GLUI_Spinner* pcSpinner = pcGlui->add_spinner_to_panel(pcPanel_Color, pszChannels[c], GLUI_SPINNER_FLOAT, &pfColor[c]);
			pcSpinner->set_float_limits(0.0f, 1.0f);
		}
	}

	pcGlui->add_checkbox_to_panel(pcPanel_Basic, "Plot Boxes?", &iIsPlottingBoxs,
				IAddWid(GLUI_EVENT_PLOT_BOXES), CGlutWin::_GluiCB_static);
	pcGlui->add_checkbox_to_panel(pcPanel_Basic, "Show Max. Prob.?", &iIsShowingMaxProb);
	pcGlui->add_checkbox_to_panel(pcPanel_Basic, "Not Recursive?", &iIsNotRecursive);

	/////////////////////////////////////////////////////////
	// setup the color editor
	cColorEditor.AddGlui(
		(CGlutWin*)this, 
		pcGlui, 
		NULL,
		pcBlockTree->uMaxLevel);

	/////////////////////////////////////////////////////////
	// setup the cluster editor
	cClusterEditor.AddGlui(
		(CGlutWin*)this, 
		pcGlui, 
		NULL,
		pcBlockTree->uMaxLevel,
		pcBlockTree->UGetNrOfBins());
}

void 
CHistView::
	_TimerFunc(unsigned short usEvent)
{
}

// ADD-BY-LEETEN 03/17/2013-BEGIN
void 
CHistView::_DisplayLevelHistograms()
{
	glTranslatef(-1.0, -1.0, 0.0f);
	glScalef(2.0, 2.0, 1.0f);
	glScalef(1.0f/(float)(i2BinRange.y + 1 - i2BinRange.x), 1.0f/(float)(iNrOfLevelsToDisplay + 1), 1.0f);
	glTranslatef(-(float)i2BinRange.x, 0.0f, 0.0f);

	////////////////////////////////////////////////////////////////
	// plot the axis
	for(size_t l = 0; l <= (size_t)iNrOfLevelsToDisplay; l++)
	{
		glPushMatrix();
		glTranslatef(0.0f, (float)l, 0.0f);

		// plot the axis
		glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
		glBegin(GL_LINES);
		glVertex2i(0, 0);
		glVertex2i((int)pcBlockTree->UGetNrOfBins(), 0);
		glEnd();

		// ADD-BY-LEETEN 02/11/2013-BEGIN
		if( iIsShowingMaxProb )
		{
			char szMaxProb[1024];
			sprintf(szMaxProb, "Prob = %f", pcBlockTree->vdLevelBinMax[l]);
			_DrawString3D(szMaxProb, (float)i2BinRange.x, 0.88f);		
		}
		// ADD-BY-LEETEN 02/11/2013-END
		glPopMatrix();
	}

	// plot the histograms
	glPushAttrib(
		GL_COLOR_BUFFER_BIT |
		0 );
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	pcBlockTree->_RenderHistograms(
		CBlock::MODE_NONE, 
		(size_t)iNrOfLevelsToDisplay,
		i2BinRange,
		this->f4DefaultColor
		);
	pcBlockTree->_RenderHistograms(
		CBlock::MODE_ASSIGNED, 
		(size_t)iNrOfLevelsToDisplay,
		i2BinRange,
		this->f4DefaultColor
		);

	// ADD-BY-LEETEN 06/23/2013-BEGIN
	for(size_t b = 0; b < CBlock::VGetBlocksRenderedByPCP().size(); b++)
	{
		CBlock::VGetBlocksRenderedByPCP()[b]->_RenderHistogram(CBlock::F4GetPCPColor(), i2BinRange, pcBlockTree->vdLevelBinMax);
	}
	// ADD-BY-LEETEN 06/23/2013-END

	glPopAttrib();
		// GL_COLOR_BUFFER_BIT |

	if( cClusterEditor.iIsActive )
	{
		// ADD-BY-LEETEN 05/07/2013-BEGIN
		glPushAttrib(
			GL_COLOR_BUFFER_BIT |
			0 );
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		pcBlockTree->_RenderHistograms(
			CBlock::MODE_SELECTED_BY_HIST,
			(size_t)iNrOfLevelsToDisplay,
			i2BinRange,
			cClusterEditor.f4Color
			);
		glPopAttrib();
			// GL_COLOR_BUFFER_BIT |
		// ADD-BY-LEETEN 05/07/2013-END
		glPushAttrib(
			GL_LINE_BIT | 
			0);
		glLineWidth(4.0f);
		glLineStipple(4, 0xCCCC);
		glEnable(GL_LINE_STIPPLE);
		glPushMatrix();
		glTranslatef(0.0f, (float)cClusterEditor.iLevel, 0.0f);
		// plot the lines
		glColor4f(cClusterEditor.f4Color.x, cClusterEditor.f4Color.y, cClusterEditor.f4Color.z, 0.1f);
		for(size_t i = 0; i < 2; i++)
		{
			glBegin(GL_LINE_STRIP);
			for(size_t b = i2BinRange.x; b <= i2BinRange.y; b++)
			{
				float fProb = (i)?this->cClusterEditor.vf2BinRanges[b].x:cClusterEditor.vf2BinRanges[b].y;
				fProb = max(min(fProb, (float)pcBlockTree->vdLevelBinMax[cClusterEditor.iLevel]), 0.0f);
				fProb /= pcBlockTree->vdLevelBinMax[cClusterEditor.iLevel];
				glVertex2f((float)b, fProb);
				glVertex2f((float)b + 1.0f, fProb);
			}
			glEnd();
		}
		// ADD-BY-LEETEN 02/03/2013-BEGIN
		// plot a rectangle to highlight the current bin
		glDisable(GL_LINE_STIPPLE);
		float l = (float)cClusterEditor.iBin;
		float r = l + 1.0f;
		float b = cClusterEditor.f2Prob.x;
		float t = cClusterEditor.f2Prob.y;
		b = max(min(b, (float)pcBlockTree->vdLevelBinMax[cClusterEditor.iLevel]), 0.0f) / pcBlockTree->vdLevelBinMax[cClusterEditor.iLevel];
		t = max(min(t, (float)pcBlockTree->vdLevelBinMax[cClusterEditor.iLevel]), 0.0f) / pcBlockTree->vdLevelBinMax[cClusterEditor.iLevel];
		glBegin(GL_LINE_STRIP);
		glVertex2f(l, b);
		glVertex2f(r, b);
		glVertex2f(r, t);
		glVertex2f(l, t);
		glVertex2f(l, b);
		glEnd();
		// ADD-BY-LEETEN 02/03/2013-END
		glPopMatrix();
		glPopAttrib(
			// GL_LINE_BIT
		);
	}
}

void 
CHistView::_DisplayFunc()
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

	_DisplayLevelHistograms();

	// restore the coordinate system back
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void 
CHistView::_KeyboardFunc(unsigned char ubKey, int x, int y)
{
}

void 
CHistView::_ReshapeFunc(int w, int h)
{
}

// ADD-BY-LEETEN 02/03/2013-BEGIN
void
CHistView::
	_GetColorEditor(
		bool& bIsValid,
		vector<size_t>& vuWaveletSub,
		vector<size_t>& vuLocalSub,
		vector<size_t>& vuGlobalSub,
		size_t& uGlobal,
		void* _Reserved
	)
{
	bIsValid = true;
	if( !cColorEditor.iIsActive )
		bIsValid = false;

	if( cColorEditor.iLevel > (int)pcBlockTree->uMaxLevel || cColorEditor.iLevel < 0 )
		bIsValid = false;

	if( !bIsValid )
		return;

	vuWaveletSub.assign(pcBlockTree->UGetNrOfDims(), cColorEditor.iLevel);

	vector<size_t> vuGlobalCoefBase, vuLocalCoefLengths;
	pcBlockTree->_ConvertWaveletSubToLevels(vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths);

	size_t uNrOfBlocks = 1;
	for(size_t d = 0; d < vuLocalCoefLengths.size(); d++)
		uNrOfBlocks *= vuLocalCoefLengths[d];

	if( cColorEditor.iBlock >= (int)uNrOfBlocks || cColorEditor.iBlock < 0 )
		bIsValid = false;

	if( !bIsValid )
		return;

	// now convert the local ID to local subscript
	WaveletSAT::_ConvertIndexToSub(cColorEditor.iBlock, vuLocalSub, vuLocalCoefLengths);

	vuGlobalSub.resize(pcBlockTree->UGetNrOfDims());
	for(size_t d = 0; d < vuLocalCoefLengths.size(); d++)
		vuGlobalSub[d] = vuGlobalCoefBase[d] + vuLocalSub[d];

	uGlobal = WaveletSAT::UConvertSubToIndex(vuGlobalSub, pcBlockTree->VGetCoefLengths());
}
// ADD-BY-LEETEN 02/03/2013-END

void 
CHistView::_GluiFunc(unsigned short usValue)
{
	switch(usValue)
	{
	case GLUI_EVENT_BIN_RANGE:
	{
		pcBlockTree->_CompMaxProb(i2BinRange);
	}	break;
	case GLUI_EVENT_PLOT_BOXES:
	{
	} break;
	case GLUI_EVENT_COLOR_ASSIGN:
	{
	} break;
	case GLUI_EVENT_COLOR_RESET:
	{
	} break;
	case GLUI_EVENT_CLUSTER_ASSIGN:
	{
		if(cClusterEditor.iIsActive)
		{
			pcBlockTree->_SelectByHist(
				CBlock::MODE_ASSIGNED, 
				true, 
				cClusterEditor.f4Color,
				i2BinRange, 
				cClusterEditor.vf2BinRanges, 
				cClusterEditor.iLevel);
		}
	} break;
	case GLUI_EVENT_CLUSTER_RESET:
	{
		if(cClusterEditor.iIsActive)
		{
			pcBlockTree->_SelectByHist(
				CBlock::MODE_NONE, 
				true, 
				cClusterEditor.f4Color,
				i2BinRange, 
				cClusterEditor.vf2BinRanges, 
				cClusterEditor.iLevel);
		}
	} break;
	case GLUI_EVENT_CLUSTER_RESET_PROB:
	{
		cClusterEditor.f2Prob = make_float2(0.0f, 1.0f);
	// MOD-BY-LEETEN 05/07/2013-FROM:	} break;
	}
	// MOD-BY-LEETEN 05/07/2013-END
	case GLUI_EVENT_CLUSTER_EDITING:
	{
		cClusterEditor.vf2BinRanges[cClusterEditor.iBin] = cClusterEditor.f2Prob;
		#if	0	// MOD-BY-LEETEN 05/07/2013-FROM:
		/*
		pcBlockTree->_ResetSelection(CBlock::MODE_SELECTED_BY_HIST);
		pcBlockTree->_SelectByHist(
			CBlock::MODE_SELECTED_BY_HIST, 
			true, 
			cClusterEditor.f4Color,
			i2BinRange, 
			cClusterEditor.vf2BinRanges, 
			cClusterEditor.iLevel);
		*/
		#else	// MOD-BY-LEETEN 05/07/2013-TO:
		if( cClusterEditor.iIsActive )
		{
			pcBlockTree->_ResetSelection(CBlock::MODE_SELECTED_BY_HIST);
			pcBlockTree->_SelectByHist(
				CBlock::MODE_SELECTED_BY_HIST, 
				true, 
				cClusterEditor.f4Color,
				i2BinRange, 
				cClusterEditor.vf2BinRanges, 
				cClusterEditor.iLevel);
		}
		#endif	// MOD-BY-LEETEN 05/07/2013-END
	} break;
	}
}

//////////////////// Constructors/Destructors //////////////////// 
CHistView::
	CHistView(void)
{
	f4DefaultColor = make_float4(0.0f, 0.0f, 0.0f, 0.1f);	// ADD-BY-LEETEN 02/11/2013

	// add a panel for the UI control
	_AddGluiSubwin(GLUI_SUBWINDOW_LEFT);
}

CHistView::
	~CHistView(void)
{
}

