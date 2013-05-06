#pragma once

#include <vector_functions.h>

#include "GlutWin.h"
#include "SepDWTHeader.h"
using namespace WaveletSAT;

#include "SATSepDWTBlockTree.h"

namespace SATSepDWT {
struct 
CPCPView:	
	virtual public CSepDWTHeader,
	virtual public CGlutWin
{
protected:
	enum {
		// event for GLUI
		GLUI_EVENT_BASE = 0x1001,
		GLUI_EVENT_FILTER_ASSIGN,
		GLUI_EVENTS_MAX_ID	
	};

	///////////////////////////////////////////////////////////////////////
	//! The current pressed button (0 means that no button is pressed now).
	float fWidth;
	float4 f4Color;

	struct CFilter 
	{
		int iTarget;
		/*
		int iLevel;
		float2 f2Prob;
		*/
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

			GLUI_Panel *pcPanel_Color = pcGlui->add_panel_to_panel(pcPanel, "Color");
			static char* pszChannels[] = {"R", "G", "B", "A"};
			float *pfColor = &f4Color.x;
			for(int c = 0; c < sizeof(pszChannels)/sizeof(pszChannels[0]); c++)
			{
				GLUI_Spinner* pcSpinner = pcGlui->add_spinner_to_panel(pcPanel_Color, pszChannels[c], GLUI_SPINNER_FLOAT, &pfColor[c]);
				pcSpinner->set_float_limits(0.0f, 1.0f);
			}
			pcGlui->add_button_to_panel(pcPanel, "Assign", pcWin->IAddWid(GLUI_EVENT_FILTER_ASSIGN), CGlutWin::_GluiCB_static);

			vf2LevelProb.assign(uNrOfLevels, make_float2(0.0f, 1.0f));
			for(size_t l = 0; l < uNrOfLevels; l++)
			{
				char szLevel[32];
				sprintf(szLevel, "Level %d", l);
				GLUI_Rollout* pcPanel_Level = pcGlui->add_rollout_to_panel(pcPanel, szLevel, 0);
				GLUI_Spinner* pcSpinner_Up = pcGlui->add_spinner_to_panel(pcPanel_Level, "Up", GLUI_SPINNER_FLOAT, &vf2LevelProb[l].y);
					pcSpinner_Up->set_float_limits(0.0f, 1.0f);
				GLUI_Spinner* pcSpinner_Low = pcGlui->add_spinner_to_panel(pcPanel_Level, "Down", GLUI_SPINNER_FLOAT, &vf2LevelProb[l].x);
					pcSpinner_Low->set_float_limits(0.0f, 1.0f);
			}
		}

		CFilter()
		{
			f4Color = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
		}
	} cFilter;

	CBlockTree *pcBlockTree;
public:

	enum {
		EVENT_ASSIGN_BLOCKS,
		NR_OF_EVENTS
	};

	void
	_SetBlockTree
	(
		CBlockTree *pcBlockTree,
		void *_Reserved = NULL
	)
	{
		this->pcBlockTree = pcBlockTree;
	}

	///////////////////////////////////////////////////////////
	// GlutWin interfaces
	void _InitFunc();
	void _GluiFunc(unsigned short usValue);
	void _DisplayFunc();

	CPCPView(void);
	~CPCPView(void);
};

}
