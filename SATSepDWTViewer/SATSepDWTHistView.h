#pragma once

#include <vector_functions.h>

#include "GlutWin.h"

#include "SATSepDWT.h"

#include "SATSepDWTPCPView.h"	// ADD-BY-LEETEN 02/14/2013

struct 
CSATSepDWTHistView:	
	virtual public CSATSepDWT,
	virtual public CGlutWin
{
protected:
	enum {
		// event for GLUI
		GLUI_EVENT_BASE = 0x901,

		GLUI_EVENT_PLOT_BOXES,	// ADD-BY-LEETEN 02/11/2013

		GLUI_EVENT_BIN_RANGE,	// ADD-BY-LEETEN 02/14/2013

		// ADD-BY-LEETEN 02/03/2013-BEGIN
		GLUI_EVENT_COLOR_ASSIGN,	
		GLUI_EVENT_COLOR_RESET,		
		// ADD-BY-LEETEN 02/03/2013-END

		// ADD-BY-LEETEN 02/03/2013-BEGIN
		GLUI_EVENT_CLUSTER_EDITING,
		GLUI_EVENT_CLUSTER_ASSIGN,	
		GLUI_EVENT_CLUSTER_RESET,		
		GLUI_EVENT_CLUSTER_RESET_PROB,
		// ADD-BY-LEETEN 02/03/2013-END

		GLUI_EVENTS_MAX_ID	
	};

	int iMaxLevel;	// The current #levels to display

	///////////////////////////////////////////////////////////////////////
	#if	0	// DEL-BY-LEETEN 02/03/2013-BEGIN
	//! The flag whether the editging is ON
	int iIsEditingCluster;	

	struct CCluster
	{
		enum{
			NR_OF_CLUSTERS_PER_LEVEL = 16,
		};

		float pfColor[4];
		vector<float2> vf2BinRanges;

		CCluster(){}

		CCluster(size_t uNrOfBins)
		{
			memset(pfColor, 0, sizeof(pfColor));
			vf2BinRanges.assign(uNrOfBins, make_float2(0.0f, 1.0f));
		}
	};
	enum 
	{
		GLUI_EVENT_CLUSTER_EDITING,
		NR_OF_GLUI_EVENTS
	};
	#endif	// DEL-BY-LEETEN 02/03/2013-END

	// ADD-BY-LEETEN 02/03/2013-BEGIN
	struct CColorEditor 
	{
		int iLevel;			// the level ID of this block 
		int iBlock;			// the block ID of this block within the current level
		int iIsActive;		// enable or disable the editing
		int iIsAssigned;	// assign the color to the current block
		float4 f4Color;

		CColorEditor():
			iLevel(0),
			iBlock(0),
			iIsActive(0)
		{
			// MOD-BY-LEETEN 02/11/2013-FROM:			memset(&f4Color, 0, sizeof(f4Color));
			f4Color = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
			// MOD-BY-LEETEN 02/11/2013-END
		}

		void
		AddGlui(
			CGlutWin* pcWin, 
			GLUI *pcGlui, 
			GLUI_Panel* pcParentPanel,
			size_t uMaxLevel,
			void* _Reserved = NULL)
		{
			GLUI_Panel *pcPanel = NULL;
			if(!pcParentPanel)
				pcPanel = pcGlui->add_rollout("Color", 1);
			else
				pcPanel = pcGlui->add_rollout_to_panel(pcParentPanel, "Color", 1);

			pcGlui->add_checkbox_to_panel(pcPanel, "Active?", &iIsActive);	// is the editor on?

			GLUI_Panel *pcPanel_Color = pcGlui->add_panel_to_panel(pcPanel, "Color");
			static char* pszChannels[] = {"R", "G", "B", "A"};
			float *pfColor = &f4Color.x;
			for(int c = 0; c < sizeof(pszChannels)/sizeof(pszChannels[0]); c++)
			{
				GLUI_Spinner* pcSpinner = pcGlui->add_spinner_to_panel(pcPanel_Color, pszChannels[c], GLUI_SPINNER_FLOAT, &pfColor[c]);
				pcSpinner->set_float_limits(0.0f, 1.0f);
			}
			GLUI_Spinner* pcSpinner_Level = pcGlui->add_spinner_to_panel(pcPanel, "Level",	GLUI_SPINNER_INT, &iLevel);	
			pcSpinner_Level->set_int_limits(0, uMaxLevel);	// 20: a big number for now

			pcGlui->add_spinner_to_panel(pcPanel, "Block",	GLUI_SPINNER_INT, &iBlock);

			pcGlui->add_button_to_panel(pcPanel, "Assign", 
				pcWin->IAddWid(GLUI_EVENT_COLOR_ASSIGN), CGlutWin::_GluiCB_static);
			pcGlui->add_button_to_panel(pcPanel, "Reset", 
				pcWin->IAddWid(GLUI_EVENT_COLOR_RESET), CGlutWin::_GluiCB_static);
		}
	} cColorEditor;
	// ADD-BY-LEETEN 02/03/2013-END

	#if	0	// MOD-BY-LEETEN 02/03/2013-FROM:
	struct CEditing {
		int iLevel;
		int iID;
		int iBin;
		CCluster cCluster;
		float2 f2Prob;

		CEditing():
			iLevel(0),
			iID(0),
			iBin(0)
		{
			f2Prob = make_float2(0.0f, 1.0f);
		}
	} cEditing;
	#else	// MOD-BY-LEETEN 02/03/2013-TO:
	struct CClusterEditor
	{
		int iIsActive;
		int iLevel;
		float4 f4Color;
		int iBin;
		float2 f2Prob;
		vector<float2> vf2BinRanges;

		void
		AddGlui(
			CGlutWin* pcWin, 
			GLUI *pcGlui, 
			GLUI_Panel* pcParentPanel,
			size_t uMaxLevel,
			size_t uNrOfBins,
			void *_Reserved = NULL)
		{
			GLUI_Panel *pcPanel = NULL;
			if(!pcParentPanel)
				pcPanel = pcGlui->add_rollout("Cluster", 0);
			else
				pcPanel = pcGlui->add_rollout_to_panel(pcParentPanel, "Cluster", 0);

			pcGlui->add_checkbox_to_panel(pcPanel, "Active?", &iIsActive);	// is the editor on?

			GLUI_Panel *pcPanel_Color = pcGlui->add_panel_to_panel(pcPanel, "Color");
			static char* pszChannels[] = {"R", "G", "B", "A"};
			float *pfColor = &f4Color.x;
			for(int c = 0; c < sizeof(pszChannels)/sizeof(pszChannels[0]); c++)
			{
				GLUI_Spinner* pcSpinner = pcGlui->add_spinner_to_panel(pcPanel_Color, pszChannels[c], GLUI_SPINNER_FLOAT, &pfColor[c]);
				pcSpinner->set_float_limits(0.0f, 1.0f);
			}

			GLUI_Panel *pcPanel_Range = pcGlui->add_panel_to_panel(pcPanel, "Range");
			GLUI_Spinner* pcSpinner_Prob;
			pcSpinner_Prob = pcGlui->add_spinner_to_panel(pcPanel_Range, "Upper",		
				GLUI_SPINNER_FLOAT, &f2Prob.y,
				pcWin->IAddWid(GLUI_EVENT_CLUSTER_EDITING), CGlutWin::_GluiCB_static);
			pcSpinner_Prob->set_float_limits(0.0f, 1.0f);
			pcSpinner_Prob = pcGlui->add_spinner_to_panel(pcPanel_Range, "Lower",		
				GLUI_SPINNER_FLOAT, &f2Prob.x,
				pcWin->	IAddWid(GLUI_EVENT_CLUSTER_EDITING), CGlutWin::_GluiCB_static);
			pcSpinner_Prob->set_float_limits(0.0f, 1.0f);

			GLUI_Spinner* pcSpinner_Level = pcGlui->add_spinner_to_panel(pcPanel, "Level",	GLUI_SPINNER_INT, &iLevel);	
			pcSpinner_Level->set_int_limits(0, uMaxLevel);

			vf2BinRanges.assign(uNrOfBins, make_float2(0.0f, 1.0f));
			GLUI_Spinner* pcSpinner_Bin = pcGlui->add_spinner_to_panel(
				pcPanel, "Bin",	GLUI_SPINNER_INT, &iBin,
				pcWin->IAddWid(GLUI_EVENT_CLUSTER_EDITING), CGlutWin::_GluiCB_static);
			// MOD-BY-LEETEN 02/14/2013-FROM:			pcSpinner_Bin->set_int_limits(0, uNrOfBins);
			pcSpinner_Bin->set_int_limits(0, uNrOfBins - 1);
			// MOD-BY-LEETEN 02/14/2013-END

			pcGlui->add_button_to_panel(pcPanel, "Assign", 
				pcWin->IAddWid(GLUI_EVENT_CLUSTER_ASSIGN), CGlutWin::_GluiCB_static);
			pcGlui->add_button_to_panel(pcPanel, "Reset", 
				pcWin->IAddWid(GLUI_EVENT_CLUSTER_RESET), CGlutWin::_GluiCB_static);
			pcGlui->add_button_to_panel(pcPanel, "Reset Prob", 
				pcWin->IAddWid(GLUI_EVENT_CLUSTER_RESET_PROB), CGlutWin::_GluiCB_static);
		}

		CClusterEditor()
		{
			f2Prob = make_float2(0.0f, 1.0f);
			f4Color = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
		}
	} cClusterEditor;
	#endif	// MOD-BY-LEETEN 02/03/2013-END

	#if	0	// DEL-BY-LEETEN 02/03/2013-BEGIN
	//! Arrays of cluster
	/*! To simplify the program, each level has its own cluster
	*/
	vector< vector<CCluster> > vvcClusters;	
	#endif	// DEL-BY-LEETEN 02/03/2013-END

	//! The current pressed button (0 means that no button is pressed now).
	int iButton;

	// MOD-BY-LEETEN 02/14/2013-FROM:	int iMinBin;
	int2 i2BinRange;
	// MOD-BY-LEETEN 02/14/2013-END

	// ADD-BY-LEETEN 02/03/2013-BEGIN
	size_t uMaxLevel;

	//! The color per coefficient whether this coefficient has its own color
	vector< pair<bool, float4> > vpairCoefColors;	
	// ADD-BY-LEETEN 02/03/2013-END

	// ADD-BY-LEETEN 02/06/2013-BEGIN
	//! A flag whether the blocks with their own color should be plotted
	int iIsPlottingBoxs;

	//! an array to store these blocks and their color
	vector< pairBlockColor > vpairBlockColors;
	// ADD-BY-LEETEN 02/06/2013-END

	// MOD-BY-LEETEN 02/14/2013-FROM:	vector< vector<WaveletSAT::typeWavelet> > vvdLevelBinMax;	// ADD-BY-LEETEN 02/06/2013
	vector<WaveletSAT::typeWavelet> vdLevelBinMax;	
	// MOD-BY-LEETEN 02/14/2013-END

	// ADD-BY-LEETEN 02/10/2013-BEGIN
	// queue of the children of its own color
	vector<size_t> vuChildrenWithColor;
	// ADD-BY-LEETEN 02/10/2013-END

	// ADD-BY-LEETEN 02/11/2013-BEGIN
	float4 f4DefaultColor;

	// if this flag is not 0, the traverse of block will stop at hte blocks of its own color
	int iIsNotRecursive;

	// if this flag is not 0, the max. prob. after iMinbin will be shown
	int iIsShowingMaxProb;
	// ADD-BY-LEETEN 02/11/2013-END

	// ADD-BY-LEETEN 02/14/2013-BEGIN
	vector< vector<double> > vvdDiffWithParent;

	void
	_ConvertToCDF
	(
		vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >& vpairBinValue,
		void *_Reserved = NULL
	);

	WaveletSAT::typeWavelet
	DCompCoefs
		(
			const vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >& vpairCDF1,
			const vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >& vpairCDF2,
			void *_Reserved = NULL
		);

	void
	CSATSepDWTHistView::
		_CompWithChild
	(
		const vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >& vpairCDF,
		const size_t uLevel, 
		const vector<size_t>& vuChildLocalSub,
		void	*_Reserved = NULL
	);

	void
	CSATSepDWTHistView::
		_CompMaxProb
	(
		void	*_Reserved = NULL
	);
	// ADD-BY-LEETEN 02/14/2013-END
public:

	CSATSepDWTPCPView cSATSepDWTPCPView;	// ADD-BY-LEETEN 02/14/2013

	// ADD-BY-LEETEN 02/06/2013-BEGIN
	enum {
		EVENT_PLOTTING_BOX,
		NR_OF_EVENTS
	};
	// ADD-BY-LEETEN 02/06/2013-END

	// ADD-BY-LEETEN 02/03/2013-BEGIN
	virtual 
	void
	_GetColorEditor(
		bool& bIsValid,
		vector<size_t>& vuWaveletSub,
		vector<size_t>& vuLocalSub,
		vector<size_t>& vuGlobalSub,
		size_t& uGlobal,
		void* _Reserved = NULL
	);

	// ADD-BY-LEETEN 02/06/2013-BEGIN
	virtual
	void
	_LoadFile
	(
		const char* szFilepath,
		void *_Reserved = NULL
	);
	// ADD-BY-LEETEN 02/06/2013-END


	virtual 
	void
	_Allocate(
		void *_Reserved = NULL
	);

	void
	_RenderHistogram
	(
		size_t uLevel,
		const vector<size_t>& vuWaveletSub,
		const vector<size_t>& vuLocalSub,
		const float4 f4Color,
		bool bIsHightLighting,
		void *_Reserved = NULL
	);
	// ADD-BY-LEETEN 02/03/2013-END

	void
	// MOD-BY-LEETEN 02/11/2013-FROM:	CSATSepDWTHistView::_RenderBlock
	_RenderBlock
	// MOD-BY-LEETEN 02/11/2013-END
	(
		size_t uLevel,
		const vector<size_t>& vuWaveletSub,
		const vector<size_t>& vuLocalSub,
		// MOD-BY-LEETEN 02/03/2013-FROM:	const float pfColor[]
		const	float4 f4Color,
		bool	bIsHightLighting,
		bool	bIsNotRecursive = false,	// ADD-BY-LEETEN 02/11/2013
		void	*_Reserved = NULL
		// MOD-BY-LEETEN 02/03/2013-END
	);

	///////////////////////////////////////////////////////////
	// GlutWin interfaces
	void _InitFunc();
	void _GluiFunc(unsigned short usValue);
	void _DisplayFunc();
	void _KeyboardFunc(unsigned char ubKey, int x, int y);
	void _ReshapeFunc(int w, int h);
	void _IdleFunc();
	void _MouseFunc(int button, int state, int x, int y);
	void _MotionFunc(int x, int y);
	void _TimerFunc(unsigned short usEvent);

	//! The method to setup advanved OpenGL features
	/*!
	This method is defined such that it can be called after _InitFunc();
	*/
	void _InitGl();

	CSATSepDWTHistView(void);
	~CSATSepDWTHistView(void);
};

