#pragma once

#include <vector_functions.h>

#include "GlutWin.h"

#include "SATSepDWT.h"

struct 
CSATSepDWTHistView:	
	virtual public CSATSepDWT,
	virtual public CGlutWin
{
protected:
	enum {
		// event for GLUI
		GLUI_EVENT_BASE = 0x901,
		GLUI_EVENTS_MAX_ID	
	};

	int iMaxLevel;	// The current #levels to display

	///////////////////////////////////////////////////////////////////////
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

	//! Arrays of cluster
	/*! To simplify the program, each level has its own cluster
	*/
	vector< vector<CCluster> > vvcClusters;	

	//! The current pressed button (0 means that no button is pressed now).
	int iButton;

	int iMinBin;
public:
	void
	CSATSepDWTHistView::_RenderBlock
	(
		size_t uLevel,
		const vector<size_t>& vuWaveletSub,
		const vector<size_t>& vuLocalSub,
		const float pfColor[]
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

