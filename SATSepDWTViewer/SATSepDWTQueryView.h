#pragma once

#include <vector_functions.h>

#include "GlutWin.h"

#include "SATSepDWT.h"

struct 
CSATSepDWTQueryView:	
	virtual public CSATSepDWT,
	virtual public CGlutWin
{
protected:
	enum {
		// event for GLUI
		GLUI_EVENT_BASE = 0xa01,

		GLUI_EVENTS_MAX_ID	
	};

	WaveletSAT::typeSum dHistogramMax;
	vector<WaveletSAT::typeSum> vdHistogram;

	///////////////////////////////////////////////////////////////////////
	enum {
		NR_OF_EVENTS
	};

public:
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

	CSATSepDWTQueryView(void);
	~CSATSepDWTQueryView(void);
};

