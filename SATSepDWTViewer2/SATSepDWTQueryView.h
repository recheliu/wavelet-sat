#pragma once

#include <vector_functions.h>

#include "GlutWin.h"
#include "SATSepDWTBlockTree.h"
#include "SATSepDWTQuery.h"

namespace SATSepDWT {
struct 
CQueryView:	
	virtual public CGlutWin
{
protected:
	enum {
		// event for GLUI
		GLUI_EVENT_BASE = 0xa01,

		GLUI_EVENT_QUERY,

		GLUI_EVENTS_MAX_ID	
	};

	WaveletSAT::typeSum dHistogramMax;
	vector<WaveletSAT::typeSum> vdHistogram;

	//! The data structure for query
	struct {
		CQuery *pcQuery;
		int iIsPrintingTiming;

		void
		_AddGlui(
			CGlutWin* pcWin, 
			GLUI *pcGlui, 
			GLUI_Panel* pcParentPanel,
			const vector<size_t>& vuDimLengths,
			void *_Reserved = NULL)
		{
			GLUI_Panel *pcPanel = NULL;
			if(!pcParentPanel)
				pcPanel = pcGlui->add_panel("Query", 0);
			else
				pcPanel = pcGlui->add_panel_to_panel(pcParentPanel, "Query", 0);

			GLUI_Panel *pcPanel_Region = pcGlui->add_panel_to_panel(pcPanel, "Region");

			static char* pszAxes[] = {"X", "Y", "Z"};
			GLUI_Panel *pcPanel_Location = pcGlui->add_panel_to_panel(pcPanel_Region, "Location");
			int *piLocation = &pcQuery->i4Location.x;
			for(int c = 0; c < min(vuDimLengths.size(), sizeof(pszAxes)/sizeof(pszAxes[0])); c++)
			{
				GLUI_Spinner* pcSpinner = 
					pcGlui->add_spinner_to_panel(pcPanel_Location, pszAxes[c], GLUI_SPINNER_INT, &piLocation[c], 
					pcWin->IAddWid(GLUI_EVENT_QUERY), CGlutWin::_GluiCB_static);
				pcSpinner->set_int_limits((int)0, (int)vuDimLengths[c] - 1);
			}
			pcGlui->add_column_to_panel(pcPanel_Region, 1);

			GLUI_Panel *pcPanel_Size = pcGlui->add_panel_to_panel(pcPanel_Region, "Size");
			int *piSize = &pcQuery->i4Size.x;
			for(int c = 0; c < min(vuDimLengths.size(), sizeof(pszAxes)/sizeof(pszAxes[0])); c++)
			{
				GLUI_Spinner* pcSpinner = 
					pcGlui->add_spinner_to_panel(pcPanel_Size, pszAxes[c], GLUI_SPINNER_INT, &piSize[c], 
					pcWin->IAddWid(GLUI_EVENT_QUERY), CGlutWin::_GluiCB_static);
				pcSpinner->set_int_limits((int)0, (int)vuDimLengths[c] - 1);
			}

			pcGlui->add_column_to_panel(pcPanel, 1);

			GLUI_Panel *pcPanel_Option = pcGlui->add_panel_to_panel(pcPanel, "Option");
			//! If this is checked, the reqion will be queried actively
			pcGlui->add_checkbox_to_panel(pcPanel_Option, "Plot?", &pcQuery->iIsPlot);	

			//! If this is checked, the query time will be printed
			pcGlui->add_checkbox_to_panel(pcPanel_Option, "Print Timing?", &iIsPrintingTiming);	

			/*
			pcGlui->add_button_to_panel(pcPanel_Option, "Query", 
				pcWin->IAddWid(GLUI_EVENT_QUERY), CGlutWin::_GluiCB_static);
			*/
		}

	} cQuery;

	CBlockTree *pcBlockTree;

public:
	///////////////////////////////////////////////////////////////////////
	enum {
		EVENT_QUERY,
		NR_OF_EVENTS
	};

	void
	_SetQuery(
		CQuery *pcQuery,
		void *_Reserved = NULL
	)
	{
		this->cQuery.pcQuery = pcQuery;
	}

	void
	_SetBlockTree
	(
		CBlockTree *pcBlockTree,
		void *_REversed = NULL
	)
	{
		this->pcBlockTree = pcBlockTree;
	}

	virtual
	void
	_SetHistgroam
	(
		const vector<WaveletSAT::typeSum>& vdHistogram,
		void *_Reversed = NULL
	)
	{
		if(	vdHistogram.size() != this->vdHistogram.size() )
			this->vdHistogram.resize(vdHistogram.size());

		copy(vdHistogram.begin(), vdHistogram.end(), this->vdHistogram.begin());

		dHistogramMax = -(WaveletSAT::typeSum)HUGE_VAL;
		for(size_t b = 0; b < vdHistogram.size(); b++)
			dHistogramMax = max(dHistogramMax, vdHistogram[b]);
	}

	virtual
	void 
	_QueryHistogram
	(
		void *_Reserved = NULL
	);

	/*
	virtual
	void
	_SetQueryLocation
	(
		const int4& i4Location,
		void *_Reserved = NULL
	);
	*/

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

	CQueryView(void);
	~CQueryView(void);
};
}