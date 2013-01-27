#include <GL/glew.h>

#include "shader.h"	
#include "libbuf3d.h"

#include "libclock.h"	

#include "SATSepDWTHistView.h"

void
CSATSepDWTHistView::_RenderBlock
(
	size_t uLevel,
	const vector<size_t>& vuWaveletSub,
	const vector<size_t>& vuLocalSub,
	const float pfColor[]
)
{
	if( (int)uLevel > iMaxLevel )
		return;

	// fetch the coordinates and plot it 
	vector<size_t> vuGlobalCoefBase;
	vector<size_t> vuLocalCoefLengths;
	this->_ConvertWaveletSubToLevels(vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths);
	vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > > vpairBinCoefs;
	_GetCoefSparse( 
		WaveletSAT::UConvertSubToIndex(vuWaveletSub,	vuDimLevels),
		WaveletSAT::UConvertSubToIndex(vuLocalSub,		vuLocalCoefLengths),
		vpairBinCoefs );
	float fSum = 0.0f;
	for(vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > >::iterator 
			ivpairBinCoef = vpairBinCoefs.begin();
		ivpairBinCoef != vpairBinCoefs.end();
		ivpairBinCoef++)
	{
		if( ivpairBinCoef->first < (size_t)iMinBin )
			continue;
		fSum += (float)ivpairBinCoef->second;
	}

	glPushMatrix();
	glTranslatef(0.0f, (float)iMaxLevel - (float)uLevel, 0.0f);
	glColor4fv(pfColor);
	glBegin(GL_LINE_STRIP);
	int iPrevBin = 0;
	for(vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > >::iterator 
			ivpairBinCoef = vpairBinCoefs.begin();
		ivpairBinCoef != vpairBinCoefs.end();
		ivpairBinCoef++)
	{
		int iBin = (int)ivpairBinCoef->first;
		if( iBin < iMinBin )
			continue;

		if( iBin - iPrevBin > 1 )
		{
			glVertex2f((float)iPrevBin + 1, 0.0f);
			glVertex2f((float)iBin, 0.0f);
		}

		float fCount = (float)ivpairBinCoef->second;
		glVertex2f((float)iBin, fCount/fSum);
		glVertex2f((float)iBin+1, fCount/fSum);
		iPrevBin = iBin;
	}
	if( iPrevBin != (int)UGetNrOfBins() - 1 )
	{
		glVertex2f((float)iPrevBin + 1, 0.0f);
		glVertex2f((float)UGetNrOfBins() - 1.0f, 0.0f);
	}
	glEnd();
	glPopMatrix();

	///////////////////////////////////////////
	// now plot the children

	// first, find the global base in the next level
	vector<size_t> vuNextWaveletSub;
	vuNextWaveletSub.resize(UGetNrOfDims());
	for(size_t d = 0; d < UGetNrOfDims(); d++)
		vuNextWaveletSub[d] = vuWaveletSub[d] + 1;
		
	// scan through each child, 
	size_t uNrOfChildren = (!uLevel)?1:((size_t)1 << UGetNrOfDims());
	vector<size_t> vuChildVol;
	vuChildVol.assign(UGetNrOfDims(), 2);
	vector<size_t> vuChildSub;
	vuChildSub.resize(UGetNrOfDims());
	for(size_t c = 0; c < uNrOfChildren; c++)
	{
		WaveletSAT::_ConvertIndexToSub(c, vuChildSub, vuChildVol);
			
		for(size_t d = 0; d < UGetNrOfDims(); d++)
			vuChildSub[d] += vuLocalSub[d] * 2;

		// if this child has its own cluster, use the color of this cluster

		// otherwise, use the passed colors
		_RenderBlock(
			uLevel + 1, 
			vuNextWaveletSub, 
			vuChildSub, 
			pfColor);
	}
}

//////////////////// CGlutWin methods //////////////////// 
void 
CSATSepDWTHistView::_IdleFunc()
{
}

void 
CSATSepDWTHistView::
	_MotionFunc
	(
		int x, 
		int y
	)
{
#if	0	
	if( GLUT_LEFT_BUTTON == iButton )
	{
		if( pcEditingBox )
		{
			pcEditingBox->ppairCorners[1].first  = ...;
			pcEditingBox->ppairCorners[1].second = ... ;

		}
	}
#endif
}

void 
CSATSepDWTHistView::_MouseFunc
	(
		int button, 
		int state, 
		int x, 
		int y
	)
{
}

void 
CSATSepDWTHistView::
	_InitGl()
{
}

void 
CSATSepDWTHistView::_InitFunc()
{
	_KeepUpdateOn();
	
	/////////////////////////////////////////////
	// set up GLUI
	GLUI *pcGlui = PCGetGluiSubwin();
	GLUI_Spinner* pcSpinner_MaxLevel = pcGlui->add_spinner("Max Level", GLUI_SPINNER_INT, &iMaxLevel);
		// decdie the max level
		size_t uMaxLevel = 0;
		for(size_t d = 0; d < UGetNrOfDims(); d++)
		{
			size_t uLevel = vuDimLevels[d] - 2;
			if( !uMaxLevel )
				uMaxLevel = uLevel;
			else
				uMaxLevel = min(uMaxLevel, uLevel);
		}
		pcSpinner_MaxLevel->set_int_limits(0, (int)uMaxLevel);
	GLUI_Spinner* pcSpinner_MinBin = pcGlui->add_spinner("Min Bin", GLUI_SPINNER_INT, &iMinBin);
		pcSpinner_MinBin->set_int_limits(0, UGetNrOfBins());

	/////////////////////////////////////////////////////////
	// setup the cluste
	GLUI_Panel *pcPanel_Cluster = pcGlui->add_panel("Cluster");
	{
		GLUI_Panel *pcPanel_Color = pcGlui->add_panel_to_panel(pcPanel_Cluster, "Color");
		static char* pszChannels[] = {"R", "G", "B", "A"};
		for(int c = 0; c < sizeof(pszChannels)/sizeof(pszChannels[0]); c++)
		{
			GLUI_Spinner* pcSpinner = pcGlui->add_spinner_to_panel(pcPanel_Color , pszChannels[c], GLUI_SPINNER_FLOAT, &cEditing.cCluster.pfColor[c]);
			pcSpinner->set_float_limits(0.0f, 1.0f);
		}
		GLUI_Spinner* pcSpinner_Level = pcGlui->add_spinner_to_panel(pcPanel_Cluster, "Level",	
				GLUI_SPINNER_INT, &cEditing.iLevel, 
				IAddWid(GLUI_EVENT_CLUSTER_EDITING), CGlutWin::_GluiCB_static);
			pcSpinner_Level->set_int_limits(0, uMaxLevel);
		GLUI_Spinner* pcSpinner_ID = pcGlui->add_spinner_to_panel(pcPanel_Cluster, "ID",		
				GLUI_SPINNER_INT, &cEditing.iID,	
				IAddWid(GLUI_EVENT_CLUSTER_EDITING), CGlutWin::_GluiCB_static);
			pcSpinner_ID->set_int_limits(0, CCluster::NR_OF_CLUSTERS_PER_LEVEL - 1);
		GLUI_Spinner* pcSpinner_Bin = pcGlui->add_spinner_to_panel(pcPanel_Cluster, "Bin",		
				GLUI_SPINNER_INT, &cEditing.iBin,
				IAddWid(GLUI_EVENT_CLUSTER_EDITING), CGlutWin::_GluiCB_static);
			pcSpinner_Bin->set_int_limits(0, UGetNrOfBins() - 1);

		GLUI_Spinner* pcSpinner_Prob;
		pcSpinner_Prob = pcGlui->add_spinner_to_panel(pcPanel_Cluster, "Upper",		
				GLUI_SPINNER_FLOAT, &cEditing.f2Prob.y,
				IAddWid(GLUI_EVENT_CLUSTER_EDITING), CGlutWin::_GluiCB_static);
			pcSpinner_Prob->set_float_limits(0.0f, 1.0f);
		pcSpinner_Prob = pcGlui->add_spinner_to_panel(pcPanel_Cluster, "Lower",		
				GLUI_SPINNER_FLOAT, &cEditing.f2Prob.x,
				IAddWid(GLUI_EVENT_CLUSTER_EDITING), CGlutWin::_GluiCB_static);
			pcSpinner_Prob->set_float_limits(0.0f, 1.0f);

		pcGlui->add_checkbox_to_panel(pcPanel_Cluster, "Editing?", &iIsEditingCluster);
	}

	cEditing.cCluster = CCluster(UGetNrOfBins());

	// allocate the clusters
	vvcClusters.resize(uMaxLevel + 1);
	for(size_t l = 0; l < vvcClusters.size(); l++)
	{
		vvcClusters[l].resize(CCluster::NR_OF_CLUSTERS_PER_LEVEL);
		for(size_t c = 0; c < vvcClusters[l].size(); c++)
			vvcClusters[l][c] = CCluster(UGetNrOfBins());
	}

}

void 
CSATSepDWTHistView::
	_TimerFunc(unsigned short usEvent)
{
}

void 
CSATSepDWTHistView::_DisplayFunc()
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

	glTranslatef(-1.0, -1.0, 0.0f);
	glScalef(2.0, 2.0, 1.0f);
	glScalef(1.0f/((float)UGetNrOfBins() - 1.0f - (float)iMinBin), 1.0f/(float)(iMaxLevel + 1), 1.0f);
	glTranslatef(-(float)iMinBin, 0.0f, 0.0f);

	////////////////////////////////////////////////////////////////
	// plot the axis
	for(size_t l = 0; l <= (size_t)iMaxLevel; l++)
	{
		glPushMatrix();
		glTranslatef(0.0f, (float)iMaxLevel - l, 0.0f);

		// plot the axis
		glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
		glBegin(GL_LINES);
		glVertex2i(0, 0);
		glVertex2i((int)UGetNrOfBins(), 0);
		glEnd();

		glPopMatrix();
	}

	////////////////////////////////////////////////////////////////// 
	// plot block histograms
	static float pfDefaultColor[4] = {0.0f, 0.0f, 0.0f, 1.0f};
	vector<size_t> vuWaveletSub;
	vuWaveletSub.assign(UGetNrOfDims(), 0);
	vector<size_t> vuLocalSub;
	vuLocalSub.assign(UGetNrOfDims(), 0);
	_RenderBlock
		(
			0, 
			vuWaveletSub,
			vuLocalSub,
			pfDefaultColor
		);

	////////////////////////////////////////////////////////////////
	// plot the editing cluster
	if( iIsEditingCluster )
	{
		glPushAttrib(
			GL_LINE_BIT | 
			0);
		glLineWidth(4.0f);
		glLineStipple(4, 0x0F0F);
		glEnable(GL_LINE_STIPPLE);
		glPushMatrix();
		glTranslatef(0.0f, (float)iMaxLevel - this->cEditing.iLevel, 0.0f);
		// plot the lines
		glColor4f(this->cEditing.cCluster.pfColor[0], this->cEditing.cCluster.pfColor[1], this->cEditing.cCluster.pfColor[2], 0.1f);
		for(size_t i = 0; i < 2; i++)
		{
			glBegin(GL_LINE_STRIP);
			for(size_t b = this->iMinBin; b < this->UGetNrOfBins(); b++)
			{
				float fProb = (i)?this->cEditing.cCluster.vf2BinRanges[b].x:this->cEditing.cCluster.vf2BinRanges[b].y;
				glVertex2f((float)b, fProb);
				glVertex2f((float)b + 1.0f, fProb);
			}
			glEnd();
		}
		glPopMatrix();
		glPopAttrib(
			// GL_LINE_BIT
		);
	}

	// restore the coordinate system back
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void 
CSATSepDWTHistView::_KeyboardFunc(unsigned char ubKey, int x, int y)
{
}

void 
CSATSepDWTHistView::_ReshapeFunc(int w, int h)
{
}

void 
CSATSepDWTHistView::_GluiFunc(unsigned short usValue)
{
	switch(usValue)
	{
	case GLUI_EVENT_CLUSTER_EDITING:
	{
		cEditing.cCluster.vf2BinRanges[cEditing.iBin] = cEditing.f2Prob;
	} break;
	}
}

//////////////////// Constructors/Destructors //////////////////// 
CSATSepDWTHistView::
	CSATSepDWTHistView(void)
{
	// add a panel for the UI control
	_AddGluiSubwin(GLUI_SUBWINDOW_LEFT);
}

CSATSepDWTHistView::
	~CSATSepDWTHistView(void)
{
}

