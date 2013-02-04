#include <GL/glew.h>

#include "shader.h"	
#include "libbuf3d.h"

#include "libclock.h"	

#include "SATSepDWTHistView.h"

// ADD-BY-LEETEN 02/03/2013-BEGIN
void
CSATSepDWTHistView::
	_Allocate(
		void *_Reserved
)
{
	CSATSepDWT::_Allocate();

	vpairCoefColors.assign(this->uNrOfCoefs, pair<bool, float4>(false, make_float4(0.0f, 0.0f, 0.0f, 1.0f)));
}

void
CSATSepDWTHistView::
	_RenderHistogram
(
	size_t uLevel,
	const vector<size_t>& vuWaveletSub,
	const vector<size_t>& vuLocalSub,
	const float4 f4Color,
	bool bIsHightLighting,
	void *_Reserved
)
{
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

	glPushAttrib(
		GL_LINE_BIT |
		0 );
	glLineWidth( (bIsHightLighting)?4.0:2.0 );
	glPushMatrix();
	// MOD-BY-LEETEN 02/03/2013-FROM:	glTranslatef(0.0f, (float)iMaxLevel - (float)uLevel, 0.0f);
	glTranslatef(0.0f, (float)uLevel, 0.0f);
	// MOD-BY-LEETEN 02/03/2013-END

	glColor4fv((float*)&f4Color);

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
	glPopAttrib();
		// GL_LINE_BIT |
}
// ADD-BY-LEETEN 02/03/2013-END

void
CSATSepDWTHistView::_RenderBlock
(
	size_t uLevel,
	const vector<size_t>& vuWaveletSub,
	const vector<size_t>& vuLocalSub,
	// MOD-BY-LEETEN 02/03/2013-FROM:	const float pfColor[]
	const float4 f4Color,
	bool	bIsHighLighting,
	void	*_Reserved
	// MOD-BY-LEETEN 02/03/2013-END
)
{
	if( (int)uLevel > iMaxLevel )
		return;

	#if	0	// MOD-BY-LEETEN 02/03/2013-FROM:
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
	#else	// MOD-BY-LEETEN 02/03/2013-TO:
	float4 f4NextColor;

	_RenderHistogram
	(
		uLevel,
		vuWaveletSub,
		vuLocalSub,
		f4Color,
		bIsHighLighting
	);
	#endif	// MOD-BY-LEETEN 02/03/2013-END
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
	
	// ADD-BY-LEETEN 02/03/2013-BEGIN
	///////////////////////////////////////////////////////////////////////
	// if this children has its own color, save it to a queue first
	// display other children that do not have its own color
	// display the children that has its own color

	// queue of the children of its own color
	vector<size_t> vuChildrenWithColor;

	vector<size_t> vuNextGlobalCoefBase, vuNextLocalCoefLengths;
	this->_ConvertWaveletSubToLevels(vuNextWaveletSub, vuNextGlobalCoefBase, vuNextLocalCoefLengths);
	// ADD-BY-LEETEN 02/03/2013-END

	for(size_t c = 0; c < uNrOfChildren; c++)
	{
		WaveletSAT::_ConvertIndexToSub(c, vuChildSub, vuChildVol);
			
		for(size_t d = 0; d < UGetNrOfDims(); d++)
			vuChildSub[d] += vuLocalSub[d] * 2;

		// ADD-BY-LEETEN 02/03/2013-BEGIN
		vector<size_t> vuGlobalSub;
		vuGlobalSub.resize(UGetNrOfDims());
		for(size_t d = 0; d < vuGlobalSub.size(); d++)
			vuGlobalSub[d] = vuNextGlobalCoefBase[d] + vuChildSub[d];
		size_t uGlobal = WaveletSAT::UConvertSubToIndex(vuGlobalSub, vuCoefLengths);

		if( this->vpairCoefColors[uGlobal].first )
		{
			vuChildrenWithColor.push_back(uGlobal);
			continue;
		}
		// ADD-BY-LEETEN 02/03/2013-END

		// otherwise, use the passed colors
		_RenderBlock(
			uLevel + 1, 
			vuNextWaveletSub, 
			vuChildSub, 
			// MOD-BY-LEETEN 02/03/2013-FROM:			pfColor);
			f4Color,
			bIsHighLighting);
			// MOD-BY-LEETEN 02/03/2013-END
	}

	// ADD-BY-LEETEN 02/03/2013-BEGIN
	vector<size_t> 
		vuBlobalCoefBase,	// not used
		vuLocalCoefLengths;	// not used
	for(size_t c = 0; c < vuChildrenWithColor.size(); c++)
	{
		this->_ConvertIndexToLevels(
			vuChildrenWithColor[c], 
			vuNextWaveletSub, 
			vuChildSub, 
			vuBlobalCoefBase,		// not used
			vuLocalCoefLengths);	// not used
			
		// if this child has its own cluster, use the color of this cluster

		// otherwise, use the passed colors
		_RenderBlock(
			uLevel + 1, 
			vuNextWaveletSub, 
			vuChildSub, 
			this->vpairCoefColors[vuChildrenWithColor[c]].second,
			false);
	}
	// ADD-BY-LEETEN 02/03/2013-END
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
		// MOD-BY-LEETEN 02/03/2013-FROM:		size_t uMaxLevel = 0;
		uMaxLevel = 0;
		// MOD-BY-LEETEN 02/03/2013-END
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

	// ADD-BY-LEETEN 02/03/2013-BEGIN
	/////////////////////////////////////////////////////////
	// setup the color editor
	cColorEditor.AddGlui(
		(CGlutWin*)this, 
		pcGlui, 
		NULL,
		uMaxLevel);
	// ADD-BY-LEETEN 02/03/2013-END

	/////////////////////////////////////////////////////////
	// setup the cluster editor
	#if	0	// MOD-BY-LEETEN 02/03/2013-FROM:
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
	#else	// MOD-BY-LEETEN 02/03/2013-TO:
	cClusterEditor.AddGlui(
		(CGlutWin*)this, 
		pcGlui, 
		NULL,
		uMaxLevel,
		UGetNrOfBins());
	#endif	// MOD-BY-LEETEN 02/03/2013-END
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
	// MOD-BY-LEETEN 02/03/2013-FROM:	static float pfDefaultColor[4] = {0.0f, 0.0f, 0.0f, 1.0f};
	float4 f4DefaultColor = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
	// MOD-BY-LEETEN 02/03/2013-END
	vector<size_t> vuWaveletSub;
	vuWaveletSub.assign(UGetNrOfDims(), 0);
	vector<size_t> vuLocalSub;
	vuLocalSub.assign(UGetNrOfDims(), 0);
	_RenderBlock
		(
			0, 
			vuWaveletSub,
			vuLocalSub,
			// MOD-BY-LEETEN 02/03/2013-FROM:	pfDefaultColor
			f4DefaultColor,
			false
			// MOD-BY-LEETEN 02/03/2013-END
		);

	// ADD-BY-LEETEN 02/03/2013-BEGIN
	if( cColorEditor.iIsActive )
	{
		bool bIsValid;
		vector<size_t> vuWaveletSub, vuLocalSub, vuGlobalSub;
		size_t uGlobal;
		_GetColorEditor(bIsValid, vuWaveletSub, vuLocalSub, vuGlobalSub, uGlobal);
		if(bIsValid)
		{
			_RenderBlock
			(
				cColorEditor.iLevel, 
				vuWaveletSub,
				vuLocalSub,
				cColorEditor.f4Color,
				true
			);
		}
	}
	// ADD-BY-LEETEN 02/03/2013-END

	////////////////////////////////////////////////////////////////
	// plot the editing cluster
	// MOD-BY-LEETEN 02/03/2013-FROM:	if( iIsEditingCluster )
	if( cClusterEditor.iIsActive )
	// MOD-BY-LEETEN 02/03/2013-END
	{
		glPushAttrib(
			GL_LINE_BIT | 
			0);
		glLineWidth(4.0f);
		glLineStipple(4, 0x0F0F);
		glEnable(GL_LINE_STIPPLE);
		glPushMatrix();
		#if	0	// MOD-BY-LEETEN 02/03/2013-FROM:
		glTranslatef(0.0f, (float)iMaxLevel - this->cEditing.iLevel, 0.0f);
		// plot the lines
		glColor4f(this->cEditing.cCluster.pfColor[0], this->cEditing.cCluster.pfColor[1], this->cEditing.cCluster.pfColor[2], 0.1f);
		#else	// MOD-BY-LEETEN 02/03/2013-TO:
		glTranslatef(0.0f, (float)cClusterEditor.iLevel, 0.0f);
		// plot the lines
		glColor4f(cClusterEditor.f4Color.x, cClusterEditor.f4Color.y, cClusterEditor.f4Color.z, 0.1f);
		#endif	// MOD-BY-LEETEN 02/03/2013-END
		for(size_t i = 0; i < 2; i++)
		{
			glBegin(GL_LINE_STRIP);
			for(size_t b = this->iMinBin; b < this->UGetNrOfBins(); b++)
			{
				// MOD-BY-LEETEN 02/03/2013-FROM:	float fProb = (i)?this->cEditing.cCluster.vf2BinRanges[b].x:this->cEditing.cCluster.vf2BinRanges[b].y;
				float fProb = (i)?this->cClusterEditor.vf2BinRanges[b].x:cClusterEditor.vf2BinRanges[b].y;
				// MOD-BY-LEETEN 02/03/2013-END
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

// ADD-BY-LEETEN 02/03/2013-BEGIN
void
CSATSepDWTHistView::
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

	if( cColorEditor.iLevel > (int)this->uMaxLevel || cColorEditor.iLevel < 0 )
		bIsValid = false;

	if( !bIsValid )
		return;

	vuWaveletSub.assign(UGetNrOfDims(), cColorEditor.iLevel);

	vector<size_t> vuGlobalCoefBase, vuLocalCoefLengths;
	this->_ConvertWaveletSubToLevels(vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths);

	size_t uNrOfBlocks = 1;
	for(size_t d = 0; d < vuLocalCoefLengths.size(); d++)
		uNrOfBlocks *= vuLocalCoefLengths[d];

	if( cColorEditor.iBlock >= (int)uNrOfBlocks || cColorEditor.iBlock < 0 )
		bIsValid = false;

	if( !bIsValid )
		return;

	// now convert the local ID to local subscript
	WaveletSAT::_ConvertIndexToSub(cColorEditor.iBlock, vuLocalSub, vuLocalCoefLengths);

	vuGlobalSub.resize(UGetNrOfDims());
	for(size_t d = 0; d < vuLocalCoefLengths.size(); d++)
		vuGlobalSub[d] = vuGlobalCoefBase[d] + vuLocalSub[d];

	uGlobal = WaveletSAT::UConvertSubToIndex(vuGlobalSub, vuCoefLengths);
}
// ADD-BY-LEETEN 02/03/2013-END

void 
CSATSepDWTHistView::_GluiFunc(unsigned short usValue)
{
	switch(usValue)
	{
	// ADD-BY-LEETEN 02/03/2013-BEGIN
	case GLUI_EVENT_COLOR_ASSIGN:
	{
		bool bIsValid;
		vector<size_t> vuWaveletSub, vuLocalSub, vuGlobalSub;
		size_t uGlobal;
		_GetColorEditor(bIsValid, vuWaveletSub, vuLocalSub, vuGlobalSub, uGlobal);
		if( !bIsValid )
			return;

		vpairCoefColors[uGlobal].first = true;
		vpairCoefColors[uGlobal].second = cColorEditor.f4Color;
	} break;

	case GLUI_EVENT_COLOR_RESET:
	{
		bool bIsValid;
		vector<size_t> vuWaveletSub, vuLocalSub, vuGlobalSub;
		size_t uGlobal;
		_GetColorEditor(bIsValid, vuWaveletSub, vuLocalSub, vuGlobalSub, uGlobal);
		if( !bIsValid )
			return;
		vpairCoefColors[uGlobal].first = false;
	} break;
	// ADD-BY-LEETEN 02/03/2013-END

	// ADD-BY-LEETEN 02/03/2013-BEGIN
	case GLUI_EVENT_CLUSTER_RESET_PROB:
		cClusterEditor.f2Prob = make_float2(0.0f, 1.0f);
		break;

	case GLUI_EVENT_CLUSTER_ASSIGN:
	case GLUI_EVENT_CLUSTER_RESET:
	{
		vector<size_t> vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths;
		vuWaveletSub.assign(UGetNrOfDims(), cClusterEditor.iLevel);
		size_t uWavelet = WaveletSAT::UConvertSubToIndex(vuWaveletSub,	vuDimLevels);
		this->_ConvertWaveletSubToLevels(vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths);
		size_t uNrOfChildren = 1;
		for(size_t d = 0; d < vuLocalCoefLengths.size(); d++)
			uNrOfChildren *= vuLocalCoefLengths[d];
		vector<size_t> vuLocalSub;
		for(size_t c = 0; c < uNrOfChildren; c++)
		{
			WaveletSAT::_ConvertIndexToSub(c, vuLocalSub, vuLocalCoefLengths);
			vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > > vpairBinCoefs;
			_GetCoefSparse( 
				uWavelet,
				c, 
				vpairBinCoefs);
			float fSum = 0.0f;
			for(size_t b = 0; b < vpairBinCoefs.size(); b++)
				fSum += vpairBinCoefs[b].second;

			vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > >::iterator ivpairBinCoef = vpairBinCoefs.begin();
			bool bIsIn = true;
			for(size_t b = this->iMinBin; b < this->UGetNrOfBins(); b++)
			{
				// update the bin
				float fProb = 0.0f;
				if( vpairBinCoefs.end() == ivpairBinCoef )
					fProb = 0.0f;
				else
				{
					if( b < ivpairBinCoef->first )
						fProb = 0.0f;
					else
					{
						fProb = ivpairBinCoef->second / fSum;
						ivpairBinCoef++;
					}
				}

				if( fProb < cClusterEditor.vf2BinRanges[b].x || fProb > cClusterEditor.vf2BinRanges[b].y )
				{
					bIsIn = false;
					break;
				}
			}
			if( !bIsIn )
				continue;

			vector<size_t> vuGlobalSub;
			vuGlobalSub.resize(UGetNrOfDims());
			for(size_t d = 0; d < vuGlobalSub.size(); d++)
				vuGlobalSub[d]  = vuGlobalCoefBase[d] + vuLocalSub[d];
			size_t uGlobal = WaveletSAT::UConvertSubToIndex(vuGlobalSub, vuCoefLengths);
			switch(usValue)
			{
				case GLUI_EVENT_CLUSTER_ASSIGN:
					this->vpairCoefColors[uGlobal].first = true;
					this->vpairCoefColors[uGlobal].second = cClusterEditor.f4Color;
					break;
				case GLUI_EVENT_CLUSTER_RESET:
					this->vpairCoefColors[uGlobal].first = false;
					break;
			}
		}
	} break;
	// ADD-BY-LEETEN 02/03/2013-END

	case GLUI_EVENT_CLUSTER_EDITING:
	{
		// MOD-BY-LEETEN 02/03/2013-FROM:		cEditing.cCluster.vf2BinRanges[cEditing.iBin] = cEditing.f2Prob;
		cClusterEditor.vf2BinRanges[cClusterEditor.iBin] = cClusterEditor.f2Prob;
		// MOD-BY-LEETEN 02/03/2013-END
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

