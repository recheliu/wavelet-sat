#include <algorithm>
#include <GL/glew.h>

#include "shader.h"	
#include "libbuf3d.h"

#include "libclock.h"	

#include "SATSepDWT.h"
#include "SATSepDWTPCPView.h"

void
CSATSepDWTPCPView::
	_BuildPolylines(
		size_t uLevel,
		size_t uLocalCoef,
		vector<size_t>& vuLocalCoefLengths,
		vector<double>& vdPolyline,
		vector<size_t>& vuPath,
		void *_Reserved
	)
{
	double dDist = this->vvdDiffWithParent[uLevel][uLocalCoef];
	vdPolyline[uLevel] = dDist;
	vuPath[uLevel] = uLocalCoef;
	if( this->vvdDiffWithParent.size() - 1 == uLevel )
	{
		vvdPolylines[uLocalCoef].assign(vdPolyline.begin(), vdPolyline.end());
		vvuPaths[uLocalCoef].assign(vuPath.begin(), vuPath.end());
		return;
	}
	vector<size_t> vuLocalSub;
	_ConvertIndexToSub(uLocalCoef, vuLocalSub, vuLocalCoefLengths);

	vector<size_t> vuChildLocalCoefLengths;
	vuChildLocalCoefLengths.assign(vuLocalCoefLengths.begin(), vuLocalCoefLengths.end());
	for(size_t d = 0; d < UGetNrOfDims(); d++)
		vuChildLocalCoefLengths[d] *= 2;
	size_t uNrOfChildren = 1 << UGetNrOfDims();
	vector<size_t> vuChildVol;
	vuChildVol.assign(UGetNrOfDims(), 2);
	vector<size_t> vuChildSub;
	for(size_t c = 0; c < uNrOfChildren; c++)
	{
		_ConvertIndexToSub(c, vuChildSub, vuChildVol);
		for(size_t d = 0; d < UGetNrOfDims(); d++)
			vuChildSub[d] += vuLocalSub[d] * 2;
		_BuildPolylines(
			uLevel + 1,
			UConvertSubToIndex(vuChildSub, vuChildLocalCoefLengths),
			vuChildLocalCoefLengths,
			vdPolyline,
			vuPath);
	}
}

void
CSATSepDWTPCPView::
	_SetDiffWithParent(
		const vector< vector<double> >& vvdDiffWithParent,
		void *_Reserved
	)
{
	size_t uNrOfLevels = vvdDiffWithParent.size();
	this->vvdDiffWithParent.resize(uNrOfLevels);
	for(size_t l = 0; l < uNrOfLevels; l++)
		this->vvdDiffWithParent[l].assign(vvdDiffWithParent[l].begin(), vvdDiffWithParent[l].end());

	// now build the line per element
	vvdPolylines.resize(this->vvdDiffWithParent[uNrOfLevels-1].size());
	vvuPaths.resize(this->vvdDiffWithParent[uNrOfLevels-1].size());

	vector<double> vdPolyline;
	vdPolyline.assign(uNrOfLevels, 0.0);
	vector<size_t> vuPath;
	vuPath.assign(uNrOfLevels, 0);
	vector<size_t> vuLocalCoefLengths;
	vuLocalCoefLengths.assign(UGetNrOfDims(), 1);
	_BuildPolylines(
		1,
		0,
		vuLocalCoefLengths,
		vdPolyline,
		vuPath
	);

	/////////////////////////////////////
}

//////////////////// CGlutWin methods //////////////////// 
void 
CSATSepDWTPCPView::_InitFunc()
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

	cFilter.f2Prob = make_float2(0.0f, 1.0f);
	cFilter.vf2LevelProb.assign(UGetNrOfLevels(), make_float2(0.0f, 1.0f));
	cFilter._AddGlui(this, pcGlui, NULL, UGetNrOfLevels());
}

void 
CSATSepDWTPCPView::_DisplayFunc()
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
	size_t uNrOfLevels = vvdDiffWithParent.size();
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
		GL_LINE_BIT |
		0 );
	glLineWidth(4.0f);
	// plot a box to highlight the current level
	glBegin(GL_LINE_LOOP);
	{
		double b =  0.0;
		double t =  1.0;
		double l =  (double)cFilter.iTarget - 0.1;
		double r =  (double)cFilter.iTarget + 0.1;
		glVertex2d(l, b);
		glVertex2d(r, b);
		glVertex2d(r, t);
		glVertex2d(l, t);
	}
	glEnd();

	glPushAttrib(
		GL_COLOR_BUFFER_BIT | 
		0 );
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);

	// plot the polylines
	glLineWidth(fWidth);
	vuFilteredPolylines.clear();
	for(size_t l = 0; l < vvdPolylines.size(); l++)
	{
		vector<double>& vdPolyline = vvdPolylines[l];
		bool bIsFiltered = true;
		for(size_t c = 1; c < vdPolyline.size(); c++)
			if( vdPolyline[c] < cFilter.vf2LevelProb[c].x || vdPolyline[c] > cFilter.vf2LevelProb[c].y )
			{
				bIsFiltered = false;
				break;
			}

		if(bIsFiltered)
		{
			vuFilteredPolylines.push_back(l);
			continue;
		}

		glBegin(GL_LINE_STRIP);
		for(size_t c = 0; c < vdPolyline.size(); c++)
		{
			glVertex2d((double)c, vdPolyline[c]);
		}
		glEnd();
	}

	glLineWidth(cFilter.fWidth);
	glColor4fv(&cFilter.f4Color.x);
	for(size_t l = 0; l < vuFilteredPolylines.size(); l++)
	{
		vector<double>& vdPolyline = vvdPolylines[vuFilteredPolylines[l]];
		glBegin(GL_LINE_STRIP);
		for(size_t c = 0; c < vdPolyline.size(); c++)
		{
			glVertex2d((double)c, vdPolyline[c]);
		}
		glEnd();
	}
	glPopAttrib();
		// GL_COLOR_BUFFER_BIT | 

	glLineWidth(4.0f);

	glBegin(GL_LINES);
	for(size_t l = 1; l < cFilter.vf2LevelProb.size(); l++)
	{
		glVertex2f((float)l, cFilter.vf2LevelProb[l].x);
		glVertex2f((float)l, cFilter.vf2LevelProb[l].y);
	}
	glEnd();
	// plot a box to highlight the current level
	glBegin(GL_LINE_LOOP);
	{
		double b =  cFilter.vf2LevelProb[cFilter.iLevel].x;
		double t =  cFilter.vf2LevelProb[cFilter.iLevel].y;
		double l =  (double)cFilter.iLevel - 0.1;
		double r =  (double)cFilter.iLevel + 0.1;
		glVertex2d(l, b);
		glVertex2d(r, b);
		glVertex2d(r, t);
		glVertex2d(l, t);
	}
	glEnd();
	glPopAttrib();
		// GL_LINE_BIT | 

	// restore the coordinate system back
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void 
CSATSepDWTPCPView::_GluiFunc(unsigned short usValue)
{
	switch(usValue)
	{
	case GLUI_EVENT_PROB:
	{
		cFilter.vf2LevelProb[cFilter.iLevel] = cFilter.f2Prob;
	} break;
	case GLUI_EVENT_RESET_PROB:
	{
		cFilter.f2Prob = make_float2(0.0f, 1.0f);
		cFilter.vf2LevelProb[cFilter.iLevel] = cFilter.f2Prob;
	} break;
	case GLUI_EVENT_FILTER_ASSIGN:
	{
		vector<size_t> vuFilterBlocks;
		for(size_t p = 0; p < vuFilteredPolylines.size(); p++)
			vuFilterBlocks.push_back(vvuPaths[vuFilteredPolylines[p]][cFilter.iTarget]);
		unique(vuFilterBlocks.begin(), vuFilterBlocks.end());
/*
		CGlutWin::_GlobalCB(
			IGetId(), 
			CGlutWin::CB_MANUAL, 
			EVENT_ASSIGN_BLOCKS, 
			cFilter.iTarget,
			(vector< size_t>*)&vuFilterBlocks,
			NULL);		
*/

		vector<size_t> vuWaveletSub, vuGlobalBase, vuLocalLengths, vuLocalSub;
		vuWaveletSub.assign(UGetNrOfDims(), cFilter.iTarget);
		_ConvertWaveletSubToLevels(vuWaveletSub, vuGlobalBase, vuLocalLengths);
		vector< pairBlockColor > vpairBlockColors;
		for(size_t b = 0; b < vuFilterBlocks.size(); b++)
		{
			size_t uLocal = vuFilterBlocks[b];
			_ConvertIndexToSub(uLocal, vuLocalSub, vuLocalLengths);

			// append the blocks with its own color to vpairBlockColors
			float4 f4Left, f4Size;
			float* pfLeft = (float*)&f4Left.x;
			float* pfSize = (float*)&f4Size.x;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				float fWaveletLength = (float)vuCoefLengths[d] / (float)vuLocalLengths[d];
				pfLeft[d] = (float)vuLocalSub[d] * fWaveletLength;
				pfSize[d] = (float)fWaveletLength;
			}
			for(size_t d = UGetNrOfDims(); d < 4; d++)
			{
				pfLeft[d] = 0.0f;
				pfSize[d] = 1.0f;
			}
			vpairBlockColors.push_back(
					make_pair<pair<float4, float4>, float4>
					(
						make_pair<float4, float4>(f4Left, f4Size), 
						cFilter.f4Color
					)
				);
		}
		CGlutWin::_GlobalCB(
			IGetId(), 
			CGlutWin::CB_MANUAL, 
			EVENT_ASSIGN_BLOCKS, 
			1,
			(vector< pairBlockColor >*)&vpairBlockColors,
			NULL);		
	} break;
	}
}

//////////////////// Constructors/Destructors //////////////////// 
CSATSepDWTPCPView::
	CSATSepDWTPCPView(void)
{
	f4Color = make_float4(0.0f, 0.0f, 0.0f, 1.0f);

	// add a panel for the UI control
	_AddGluiSubwin(GLUI_SUBWINDOW_LEFT);
}

CSATSepDWTPCPView::
	~CSATSepDWTPCPView(void)
{
}

