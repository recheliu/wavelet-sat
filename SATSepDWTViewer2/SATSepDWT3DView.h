#pragma once

#include <vector_functions.h>

#include "libdvrsuite/DvrSuiteWin.h"

#include "SATSepDWTBlockTree.h"
#include "SATSepDWTQuery.h"

namespace SATSepDWT {
struct 
C3DView:	
	virtual public CDvrSuiteWin
	/*
	virtual public CClipVolume
	*/
{
protected:
	enum {
		// event for GLUI
		GLUI_EVENT_BASE = 0x901,
		GLUI_EVENTS_MAX_ID	
	};

	// The modelview matrix after the scaling
	TMatrix tModifiedModelviewMatrix;

	struct {
		bool bActive;
		double pdViewRayStep_obj[3];
		double pdCoord_obj[3];
	} cCursor3D;

	CBlockTree *pcBlockTree;
	CQuery *pcQuery;
public:
	enum {
		// event for GLUI
		EVENT_BASE = 0x900,
		EVENT_CURSOR_3D,
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

	void
	_SetQuery(
		CQuery *pcQuery,
		void *_Reserved = NULL
	)
	{
		this->pcQuery = pcQuery;
	}


	void _InitFunc();
	void _GluiFunc(unsigned short usValue);

	void _DisplayFunc();
	void _KeyboardFunc(unsigned char ubKey, int x, int y);
	void _ReshapeFunc(int w, int h);
	void _IdleFunc();
	void _MouseFunc(int button, int state, int x, int y);
	void _TimerFunc(unsigned short usEvent);
	////////// volume rendering methods //////////////
	void _BeginDisplay();
	void _EndDisplay();
	void _RenderSlab(
		int iSlab, int iNrOfSlabs,
		double pdModelviewMatrix[], double pdProjectionMatrix[], int piViewport[],
		double dMinX, double dMaxX, 
		double dMinY, double dMaxY, 
		double dMinZ, double dMaxZ);
	//! The method to setup advanved OpenGL features
	/*!
	This method is defined such that it can be called after _InitFunc();
	*/
	void _InitGl();

	C3DView(void);
	~C3DView(void);
};
}

/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
