#pragma once

#include <vector_functions.h>

#include "GlutWin.h"
#include "SepDWTHeader.h"
using namespace WaveletSAT;

struct 
CSATSepDWTPCPView:	
	virtual public CSepDWTHeader,
	virtual public CGlutWin
{
protected:
	enum {
		// event for GLUI
		GLUI_EVENT_BASE = 0x1001,
		GLUI_EVENT_PROB,
		GLUI_EVENT_RESET_PROB,
		GLUI_EVENT_FILTER_ASSIGN,
		GLUI_EVENTS_MAX_ID	
	};

	///////////////////////////////////////////////////////////////////////
	//! The current pressed button (0 means that no button is pressed now).
	vector< vector<double> > vvdDiffWithParent;

	vector< vector<double> > vvdPolylines;

	vector< vector<size_t> > vvuPaths;

	vector<size_t> vuFilteredPolylines;

	float fWidth;
	float4 f4Color;

	struct CFilter 
	{
		int iTarget;
		int iLevel;
		float2 f2Prob;
		vector< float2 > vf2LevelProb;
		float4 f4Color;
		float fWidth;

		void
		_AddGlui(
			CGlutWin* pcWin, 
			GLUI *pcGlui, 
			GLUI_Panel* pcParentPanel,
			size_t uNrOfLevels,
			void* _Reserved = NULL)
		{
			GLUI_Panel* pcPanel = pcGlui->add_rollout("Filter");
			GLUI_Spinner* pcSpinner_Target = pcGlui->add_spinner_to_panel(pcPanel, "Target", GLUI_SPINNER_INT, &iTarget);
				pcSpinner_Target->set_int_limits(1, uNrOfLevels - 1);
			GLUI_Spinner* pcSpinner_Level = pcGlui->add_spinner_to_panel(pcPanel, "Level", GLUI_SPINNER_INT, &iLevel, pcWin->IAddWid(GLUI_EVENT_PROB), CGlutWin::_GluiCB_static);
				pcSpinner_Level->set_int_limits(1, uNrOfLevels - 1);
			GLUI_Spinner* pcSpinner_Up = pcGlui->add_spinner_to_panel(pcPanel, "Up", GLUI_SPINNER_FLOAT, &f2Prob.y, pcWin->IAddWid(GLUI_EVENT_PROB), CGlutWin::_GluiCB_static);
				pcSpinner_Up->set_float_limits(0.0f, 1.0f);
			GLUI_Spinner* pcSpinner_Low = pcGlui->add_spinner_to_panel(pcPanel, "Low", GLUI_SPINNER_FLOAT, &f2Prob.x, pcWin->IAddWid(GLUI_EVENT_PROB), CGlutWin::_GluiCB_static);
				pcSpinner_Low->set_float_limits(0.0f, 1.0f);
			
			pcGlui->add_button_to_panel(pcPanel, "Reset Prob", pcWin->IAddWid(GLUI_EVENT_RESET_PROB), CGlutWin::_GluiCB_static);

			GLUI_Panel *pcPanel_Color = pcGlui->add_panel_to_panel(pcPanel, "Color");
			static char* pszChannels[] = {"R", "G", "B", "A"};
			float *pfColor = &f4Color.x;
			for(int c = 0; c < sizeof(pszChannels)/sizeof(pszChannels[0]); c++)
			{
				GLUI_Spinner* pcSpinner = pcGlui->add_spinner_to_panel(pcPanel_Color, pszChannels[c], GLUI_SPINNER_FLOAT, &pfColor[c]);
				pcSpinner->set_float_limits(0.0f, 1.0f);
			}
			GLUI_Spinner* pcSpinner_Width = pcGlui->add_spinner_to_panel(pcPanel, "Width", GLUI_SPINNER_FLOAT, &fWidth);
				pcSpinner_Width->set_float_limits(1.0f, 16.0f);

			pcGlui->add_button_to_panel(pcPanel, "Assign", pcWin->IAddWid(GLUI_EVENT_FILTER_ASSIGN), CGlutWin::_GluiCB_static);
		}

		CFilter()
		{
			f4Color = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
		}
	} cFilter;

	size_t UGetNrOfLevels()
	{
		return vvdDiffWithParent.size();
	}
public:

	enum {
		EVENT_ASSIGN_BLOCKS,
		NR_OF_EVENTS
	};

	void
	_BuildPolylines(
		size_t uLevel,
		size_t uLocalCoef,
		vector<size_t>& vuLocalCoefLengths,
		vector<double>& vdPolyline,
		vector<size_t>& vuPath,
		void *_Reserved = NULL
	);

	void
	_SetDiffWithParent(
		const vector< vector<double> >& vvdDiffWithParent,
		void *_Reserved = NULL
	);

	///////////////////////////////////////////////////////////
	// GlutWin interfaces
	void _InitFunc();
	void _GluiFunc(unsigned short usValue);
	void _DisplayFunc();

	CSATSepDWTPCPView(void);
	~CSATSepDWTPCPView(void);
};

