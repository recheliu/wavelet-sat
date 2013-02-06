#include <GL/glew.h>

#include "shader.h"	
#include "libbuf3d.h"

#include "libclock.h"	

#include "SATSepDWT3DView.h"

// ADD-BY-LEETEN 02/05/2013-BEGIN
template<typename T>
void
CSATSepDWT3DView::
	_ConvertDataToTexture
	(
		Nrrd *nin
	)
{
	T *data = (T*)nin->data;	

	size_t uNrOfDims = nin->dim;
	size_t uNrOfElements = 1;
	for(size_t d = 0; d < uNrOfDims; d++)
		uNrOfElements *= nin->axis[d].size;

	// search for the range
	dValueMin = HUGE_VAL;
	dValueMax = -HUGE_VAL;
	for(size_t v = 0; v < uNrOfElements; v++)
	{
		dValueMin = min(dValueMin, (double)data[v]);
		dValueMax = max(dValueMax, (double)data[v]);
	}
	LOG_VAR(dValueMin);
	LOG_VAR(dValueMax);

	size_t uNrOfTfEntries = vfHist.size();
	for(size_t v = 0; v < uNrOfElements; v++)
	{
		size_t uTfEntry = (size_t)floor((double)(uNrOfTfEntries - 1) * ((double)data[v] - dValueMin) / (double)(dValueMax - dValueMin));
		uTfEntry = min(max(uTfEntry, 0), uNrOfTfEntries - 1);	
		vfHist[uTfEntry] += 1.0f;
	}

	// normalize the histogram
	fMaxCount = -(float)HUGE_VAL;
	for(size_t b = 0; b < vfHist.size(); b++) 
		fMaxCount = max(fMaxCount, vfHist[b]);

	for(size_t b = 0; b < vfHist.size(); b++) 
		vfHist[b] /= fMaxCount;

	// normalize the volume to [0, 1] for floating point data
	if( nrrdTypeFloat == nin->type )
	{
		float *data = (float*)nin->data;	
		for(size_t v = 0; v < uNrOfElements; v++)
			data[v] = (data[v] - (float)dValueMin)/(float)(dValueMax - dValueMin);
	}
}

void
CSATSepDWT3DView::
	_LoadData
	(
		char* szFilepath,
		void* _Reserved
	)
{
	nin = nrrdNew();
	if (nrrdLoad(nin, szFilepath, NULL)) {
		char *err = biffGetDone(NRRD);
		LOG_ERROR(fprintf(stderr, "%s", err));
		free(err);
		exit(EXIT_FAILURE);
	}

	size_t uNrOfTfEntries = this->pcSATSepDWT->UGetNrOfBins();
	vfHist.assign(uNrOfTfEntries, 0.0f);

	switch(nin->type)
	{
	case nrrdTypeUChar:	_ConvertDataToTexture<unsigned char>(nin);	break;
	case nrrdTypeChar:	_ConvertDataToTexture<char>(nin);			break;
	case nrrdTypeShort:	_ConvertDataToTexture<short>(nin);			break;
	case nrrdTypeUShort:_ConvertDataToTexture<unsigned short>(nin);	break;
	case nrrdTypeInt:	_ConvertDataToTexture<int>(nin);			break;
	case nrrdTypeUInt:	_ConvertDataToTexture<unsigned int>(nin);	break;
	case nrrdTypeFloat:	_ConvertDataToTexture<float>(nin);			break;
	default:
		break;
	}

	// setup the transfer function
	vf4TransFunc.assign(uNrOfTfEntries, make_float4(0.0f, 0.0f, 0.0f, 0.0f));

	cTransFunc._LoadRainBow();
	cTransFunc._ExportColorMap((float*)&vf4TransFunc.front(), vf4TransFunc.size());
	cTransFunc._SetTfDomain((float)dValueMin, (float)dValueMax);
}
// ADD-BY-LEETEN 02/05/2013-END

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

// ADD-BY-LEETEN 02/05/2013-BEGIN
void 
CSATSepDWT3DView::
	_InitTf()
{
	///////////////////////////////////////////////////////////////
	// create the TF window beforehand 
	cTfWin._SetTransFunc(&cTransFunc);
	cTfWin._SetNrOfEntries(256);
	cTfWin.ICreate("Transfer Function");
	cTfWin._Set();
	cTfWin._SetHistogram(
		&vfHist.front(), 
		vfHist.size(), 
		fMaxCount, 
		dValueMin, 
		dValueMax);
	cTfWin._KeepUpdateOn();

	cTfUi._SetTransFunc(&cTransFunc);
	cTfUi.ICreate("Transfer Function Editor");
	cTfUi._SetHistogramAsBackground(
		&vfHist.front(), 
		vfHist.size(), 
		dValueMin, 
		dValueMax);
}
// ADD-BY-LEETEN 02/05/2013-END

void 
CSATSepDWT3DView::_InitFunc()
{
	// ADD-BY-LEETEN 02/05/2013-BEGIN
	_InitTf();
	_Set();	// reset the current window to this window
	// ADD-BY-LEETEN 02/05/2013-END

	CDvrWin2::_InitFunc();

	_KeepUpdateOn();

	// ADD-BY-LEETEN 02/05/2013-BEGIN
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	_DisableVerticalSync();
	_DisplayFpsOn();		// the FPS will be displayed 

	pidRayIntegral = CSetShadersByString(
		NULL
		,
		#include "ray_integral.frag.h"	
	);
	assert( pidRayIntegral );	

	////////////////////////////////////////////////////////////////
	CClipVolume::_InitFunc();

	////////////////////////////////////////////////////////////////
	// setup GLUI 
	GLUI *pcGlui = PCGetGluiWin();

	GLUI_Spinner *pcSpinner_NrOfSlices = PCGetGluiWin()->add_spinner("#Slices", GLUI_SPINNER_INT, &iNrOfSlices);	
		pcSpinner_NrOfSlices->set_int_limits(1, 4096);

						// create a spinner to control the brightness gain 
	GLUI_Spinner *pcSpinner_ThicknessGain = PCGetGluiWin()->add_spinner("Thickness Gain", GLUI_SPINNER_FLOAT, &fThicknessGain);	
	pcSpinner_ThicknessGain->set_float_limits(0.0f, 4096.0f);

	////////////////////////////////////////////////////////////////
	// upload the volume as a texture
	int iType;
	switch(nin->type)
	{
	case nrrdTypeUChar:	iType = GL_UNSIGNED_BYTE;	break;
	case nrrdTypeChar:	iType = GL_BYTE;			break;
	case nrrdTypeShort:	iType = GL_SHORT;			break;
	case nrrdTypeUShort:iType = GL_UNSIGNED_SHORT;	break;
	case nrrdTypeInt:	iType = GL_INT;				break;
	case nrrdTypeUInt:	iType = GL_UNSIGNED_INT;	break;
	case nrrdTypeFloat:	iType = GL_FLOAT;			break;
	}

	// upload the data as a texture
	switch(nin->dim)
	{
	case 3:
		_SetVolume(
			GL_LUMINANCE32F_ARB, 
			nin->data, 
			iType, 
			GL_LUMINANCE, 
			(int)nin->axis[0].size, 
			(int)nin->axis[1].size, 
			(int)nin->axis[2].size);
		break;
	case 2:
		// not supported yet. It requires the overloading of _SetVolume().		
		break;
	}

	////////////////////////////////////////////////////////////////
	_SetTransferFunc((float*)&vf4TransFunc.front(), GL_RGBA, GL_FLOAT, vf4TransFunc.size());
	_LoadSavedMatrix();		
	_SetDataValue((float)dValueMin, (float)dValueMax);
	_SetTfDomain((float)dValueMin, (float)dValueMax);
	// ADD-BY-LEETEN 02/05/2013-END
}

void 
CSATSepDWT3DView::
	_TimerFunc(unsigned short usEvent)
{
}

void 
CSATSepDWT3DView::_DisplayFunc()
{
	#if	0	// MOD-BY-LEETEN 02/05/2013-FROM:
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	#else	// MOD-BY-LEETEN 02/05/2013-TO:
	CDvrWin2::_DisplayFunc();
	#endif	// MOD-BY-LEETEN 02/05/2013-END
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
	CClipVolume::_ReshapeFunc(w, h);	// ADD-BY-LEETEN 02/05/2013
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
	// ADD-BY-LEETEN 02/05/2013-BEGIN
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	glPushMatrix();
	float fMaxDim = max(iXDim, max(iYDim, iZDim));
	glScalef(
		iXDim / fMaxDim,
		iYDim / fMaxDim,
		iZDim / fMaxDim);

	glutWireCube(2.0);	// plot the bounding box

	CClipVolume::_Create();

	glUseProgramObjectARB(	pidRayIntegral);
	SET_1F_VALUE_BY_NAME(	pidRayIntegral, "fWindowWidth",				(float)piViewport[2]);
	SET_1F_VALUE_BY_NAME(	pidRayIntegral, "fWindowHeight",			(float)piViewport[3]);
	SET_1F_VALUE_BY_NAME(	pidRayIntegral, "fThicknessGain",			fThicknessGain);

	SET_1I_VALUE_BY_NAME(	pidRayIntegral, "t2dPrevLayer",				0);
	SET_1I_VALUE_BY_NAME(	pidRayIntegral, "t3dVolume",				1);
	SET_1I_VALUE_BY_NAME(	pidRayIntegral, "t1dTf",					2);
	SET_1I_VALUE_BY_NAME(	pidRayIntegral, "t2dClipVolume",			4);
	SET_1I_VALUE_BY_NAME(	pidRayIntegral, "t2dsDepth",				5);

	SET_1F_VALUE_BY_NAME(	pidRayIntegral, "fTfDomainMin",				fTfDomainMin);
	SET_1F_VALUE_BY_NAME(	pidRayIntegral, "fTfDomainMax",				fTfDomainMax);
	SET_1F_VALUE_BY_NAME(	pidRayIntegral, "fDataValueMin",			fDataValueMin);
	SET_1F_VALUE_BY_NAME(	pidRayIntegral, "fDataValueMax",			fDataValueMax);

	glUseProgramObjectARB(0);

						// bind the volume, range, and the lookup table as textures
	glActiveTexture(GL_TEXTURE0 + 1);
	glBindTexture(GL_TEXTURE_3D, t3dVol);

	glActiveTexture(GL_TEXTURE0 + 2);
	glBindTexture(GL_TEXTURE_1D, t1dTf);

	glActiveTexture(GL_TEXTURE0 + 3);
	glBindTexture(CDvrWin2::cColor.eTarget, CDvrWin2::cColor.t2d);

	glActiveTexture(GL_TEXTURE0 + 4);
	glBindTexture(CClipVolume::cTexture.eTarget, CClipVolume::cTexture.t2d);

	glActiveTexture(GL_TEXTURE0 + 5);
	glBindTexture(CDvrWin2::cDepth.eTarget, CDvrWin2::cDepth.t2d);

	glActiveTexture(GL_TEXTURE0);
	// ADD-BY-LEETEN 02/05/2013-END
}

void 
CSATSepDWT3DView::
	_EndDisplay()
{
	// ADD-BY-LEETEN 02/05/2013-BEGIN
	glUseProgramObjectARB(0);
	glPopMatrix();
	// ADD-BY-LEETEN 02/05/2013-END
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
	// ADD-BY-LEETEN 02/05/2013-BEGIN
	glPushAttrib(
		GL_ENABLE_BIT |
		GL_DEPTH_BUFFER_BIT |
		0 );
	glEnable( GL_BLEND );
	glEnable(GL_DEPTH_TEST);

	glUseProgramObjectARB(pidRayIntegral);
	CDvrWin2::_RenderSlab(
		iSlab, iNrOfSlabs, 

		pdModelviewMatrix, pdProjectionMatrix, piViewport,
		
		dMinX, dMaxX, 
		dMinY, dMaxY, 
		dMinZ, dMaxZ);

	glPopAttrib();	// glPushAttrib(GL_DEPTH_BUFFER_BIT);
	// ADD-BY-LEETEN 02/05/2013-END
}

//////////////////// Constructors/Destructors //////////////////// 
CSATSepDWT3DView::
	CSATSepDWT3DView(void)
{
	// ADD-BY-LEETEN 02/05/2013-BEGIN
	ibIsFboEnabled = 1;							// enable the rendering to FBO

	_SetInternalColorFormat(GL_RGBA32F_ARB);	// set the depths of each chanel of the FBO as 32 bits 
	_SetInternalDepthFormat(GL_DEPTH_COMPONENT);

	fThicknessGain = 1.0f;
	iNrOfSlices = 128;

	_AddGluiWin();	// with a separate GLUI window
	// ADD-BY-LEETEN 02/05/2013-END
}

CSATSepDWT3DView::
	~CSATSepDWT3DView(void)
{
}
