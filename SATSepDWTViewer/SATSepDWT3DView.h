#pragma once

#include <vector_functions.h>

#if	0	// MOD-BY-LEETEN 02/10/2013-FROM:
#include "libclip/ClipVolume.h"
#include "libdvr2/DvrWin2.h"
// ADD-BY-LEETEN 02/05/2013-BEGIN
#include "libtfw/TfWin.h"
#include "libtfw/TfUi.h"	
// ADD-BY-LEETEN 02/05/2013-END
#include "libtfw/TransFunc.h"
#else	// MOD-BY-LEETEN 02/10/2013-TO:
#include "libdvrsuite/DvrSuiteWin.h"
#endif	// MOD-BY-LEETEN 02/10/2013-END

#include "SATSepDWT.h"

struct 
CSATSepDWT3DView:	
	#if	0	// MOD-BY-LEETEN 02/10/2013-FROM:
	virtual	public CTfUi::CReceiver,	// ADD-BY-LEETEN 02/06/2013
	virtual public CClipVolume,	// ADD-BY-LEETEN 02/05/2013
	virtual public CDvrWin2
	#else	// MOD-BY-LEETEN 02/10/2013-TO:
	virtual public CDvrSuiteWin
	#endif	// MOD-BY-LEETEN 02/10/2013-END
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

	CSATSepDWT *pcSATSepDWT;

	#if	0	// DEL-BY-LEETEN 02/10/2013-BEGIN
	// ADD-BY-LEETEN 02/05/2013-BEGIN
	//////////////////////////////////////////////////////
	// fields for volume rendering
	Nrrd *nin;

	double dValueMin;
	double dValueMax;

	// histogram
	vector<float> vfHist;
	float fMaxCount; // the maximal count of the histogram	

	GLuint pidRayIntegral;
	// ADD-BY-LEETEN 02/05/2013-END
	#endif	// DEL-BY-LEETEN 02/10/2013-END

	vector< pairBlockColor > vpairBlockColors;	// ADD-BY-LEETEN 02/06/2013

public:
	// ADD-BY-LEETEN 02/05/2013-BEGIN
	#if	0	// DEL-BY-LEETEN 02/10/2013-BEGIN
	vector<float4> vf4TransFunc;
	CTransFunc cTransFunc;

	CTfWin	cTfWin;
	CTfUi	cTfUi;

	// ADD-BY-LEETEN 02/06/2013-BEGIN
	virtual 
	void
	_Receive(
		CTransFunc *pcTransFunc
	);
	// ADD-BY-LEETEN 02/06/2013-END
	#endif	// DEL-BY-LEETEN 02/10/2013-END

	// ADD-BY-LEETEN 02/06/2013-BEGIN
	void
	_SetBlockColors
	(
		bool bIsPlottingBlocks,
		const vector< pairBlockColor >& vpairBlockColors,
		void* _Reserved = NULL
	);
	// ADD-BY-LEETEN 02/06/2013-END

	#if	0	// DEL-BY-LEETEN 02/10/2013-BEGIN
	void 
	_InitTf
	(
	);
	#endif	// DEL-BY-LEETEN 02/10/2013-END
	// ADD-BY-LEETEN 02/05/2013-END

	void
	_SetData
	(
		CSATSepDWT *pcSATSepDWT
	)
	{
		this->pcSATSepDWT = pcSATSepDWT;
	}

	#if	0	// DEL-BY-LEETEN 02/10/2013-BEGIN
	// ADD-BY-LEETEN 02/05/2013-BEGIN
	template<typename T>
	void
	_ConvertDataToTexture
	(
		Nrrd *nin
	);

	void
	_LoadData
	(
		char* szFilepath,
		void* _Reserved = NULL
	);
	// ADD-BY-LEETEN 02/05/2013-END
	#endif	// DEL-BY-LEETEN 02/10/2013-END
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

	CSATSepDWT3DView(void);
	~CSATSepDWT3DView(void);
};

/*

$Log: not supported by cvs2svn $
Revision 1.18  2012-03-30 05:54:16  leeten

[03/29/2012]
<SOM.cpp>
<SOM.h>
1. [ADD] Define methods _CompBasicComplexity() and _CompComplexity() to compute the complexity of all cluster in batch, which can be accelerated by GPUs.
2. [DEL] Skip the computation of time difference volume.

<SOMTree.cpp>
1. [ADD] Include POMCuda.h in order to use CUDA to compute the complexity.
2. [ADD] Print the timing of CSomTree::_SortClusters() when the preprocessor PRINT_CSomTree_SortClusters_TIMING is 1.

<SOMTreeView.cpp>
<SOMTreeView.h>
1. [ADD] Compute the POMs of all clusters on CPU. A new button is added to triger the event. The ID is GLUI_EVENT_COMPUTE_ALL_SOMS.

<SOMViewer_main.cpp>
1. [MOD] Change the window shape of the embedding view to squared.

Revision 1.17  2012-03-26 00:21:13  leeten

[03/25/2012]
<SOM.cpp>
<SOM.h>
1. [ADD] Add a new metric called divergence (Div).
2. [MOD] Correct the name COMPLEXOITY_* to COMPLEXITY_*.
3. [MOD] Move the computation for coordinate time series to PrecomputeStatistic(). Also remove the variable p3Df4BasisDiff
4. [ADD] Adda new entry DataName to the .ini file. If DataName is Isabel, the background textrue is loaded.

<SOMTree.cpp>
1. [MOD] Correct the name COMPLEXOITY_* to COMPLEXITY_*.

<SOMTreeNode.cpp>
1. Plot the nodes as white bounding boxes.

<SOMTreeView.cpp>
<SOMTreeView.h>
1. [ADD] Define a new method _PresetInitQuadtreeThreshold() to specify the initial qaudtree threshold.
2. [ADD] Plot the background for Isabel.
3. [MOD] Plot the quadtree nodes as black cells.

<SOMViewer_main.cpp>
1. [ADD] Add a new argument --init-quadtree-threshold tp specify the inital quadtree threshold.

Revision 1.16  2012-03-23 23:54:49  leeten

[03/23/2012]
<SOM.cpp>
<SOM.h>
1. [ADD] The MATLAB engine is enable only if this processor USE_MATLAB_FOR_EMBEDDING is 1.

<SOMEmbeddingView.cpp>
1. [ADD] The MATLAB engine is enable only if this processor USE_MATLAB_FOR_EMBEDDING is 1.
2. [MOD] Make the the occupy test more precise. The each tile in the domain stores a list of the occupied POMs, and thus the new one can precisely check whether it is intersect with the plotted one.
3. [ADD] Plot the top but occupied POM as blue bounding boxes.

<SOMTree.cpp>
1. [ADD] Implement the Telta & van Wijlk's algorithm buy using a hash table for the edges. The algorithm is enabled if the preprocessor HASH_EDGES is 1.
2. [ADD] Deump the performance of hierarchical cluster if the preprocessor PRINT_Hierarhical_Clustering_Timing is 1.
3. [ADD] Use MATLAB to do embedding if USE_MATLAB_FOR_EMBEDDING is 1.

<SOMTreeView.cpp>
<SOMTreeView.h>
1. [MOD] Set the inital quadtree threshold as 0.1.
2. [ADD] Plot the domain in Space-Only View with the data's aspect ratio. The domain can be scaled by a GLUI spiiner.
3. [MOD] Change the background color to white.
4. [ADD] Adda button to reset the clicked cluster to NULL.

<SOMViewer_main.cpp>
1. [ADD] initialize the quadtree with threshold 0.1

<SOMTree.h>
1. [ADD] Use MATLAB to do embedding if USE_MATLAB_FOR_EMBEDDING is 1.

Revision 1.15  2012-03-19 02:57:04  leeten

[03/18/2012]
<SOM.h>
<SOMBasisView.h>
<SOMCluster.h>
<SOMEmbeddingView.h>
<SOMTreeView.h>
1. [DEL] Remove old deleted code segments and labels before 03/20/2012.

Revision 1.14  2012-03-19 01:37:23  leeten

[03/18/2012]
<SOM.h>
<SOM.cpp>
1. [ADD] Add a method FGetBasicNorm() to compute the basic cluster complexity.
2. [ADD] Add a method _PrecomputeStatistics() to compute the metrics needed for fast computation of cluster complexity.
3. [MOD] In the method FGetDiffNorm(), use the FGetBasicNorm() to compute the cluster complexity and then form coumpund complexity when specfied.
4. [ADD] Read the DCT responses in _ReadSom().

<SOMBasisView.cpp>
1. [MOD] Set iIsUsingPreComputedMap as 1 so the CPU-based rendering is used by default.

<SOMEmbeddingView.cpp>
1. [MOD] Use white background instead. To make it work, the dst. factor for blending is changed to GL_ONE_MINUS_SRC_ALPHA.
2. [MOD] Use the min. of #rendering POM and the number of clusters as the threshold.
3. [MOD] The #rendering POM is used as a threshold of the rank.
4. [MOD] The prioritied selection can work for both overlapping and non-overlapping modes.
5. [MOD] Enable the correct texture target.
6. [MOD] Change the point size to 2.0.
7. [MOD] Change the color of the highlight box from white to red.

<SOMPcpView.cpp>
1. [ADD] Plot the axis of value 0.

<SOMTree.cpp>
<SOMTree.h>
1. [ADD] Add a method _SortClusters() to sort the clusters by their complexity and decide the range.
2. [MOD] Change the variable fMax1stOrderDiff to two variables fMinClusterComplexity and fMaxClusterComplexity.
3. [ADD] Add enum for the complexity metrics, and GLUI controls to specify the complexity metrics.

<SOMTreeCluster.cpp>
1. [MOD] Decide the color based on both end of the complexity range.

<SOMTreeView.h>
<SOMTreeView.cpp>
1. [ADD] Use both ends of the complexity range to generate the color.
2. [ADD] Add GLUI control to indicate the change of complexity metric.

<SOMViewer_main.cpp>
1. [MOD] Change the window size of BasisView and EmbeddingView.

Revision 1.13  2012-03-17 19:11:29  leeten

[03/17/2012]
<SOM.cpp>
1. [ADD] Compute the matrix for fast computation of 1st order and 2nd order difference.
2. [DEL] Remvoe the code to dump the time differnce field.
3. [ADD] Compute the average time difference from all grids.
4. [ADD[ Compute the histogram of time difference.

<SOMBasisView.cpp>
1. [ADD] Move the POM to the center of the windows.
2. [ADD] Remove the code to plot the time difference.

<SOMCluster.cpp>
1. [ADD] Compute the 1st order difference, 2nd order difference and its variance for the clickeding cluster.

<SOMPcpView.h>
<SOMPcpView.cpp>
1. [ADD] Add GLUI control to decide whether all constraints are shown and wherther all polylines are shown.

<SOMTree.h>
<SOMTree.cpp>
1. [ADD] Use an 1D array to store the sorted clusters with qsort. (The reason is that the stdLLlist::sort() crashed in Debug mode.

<SOMTreeView.h>
<SOMTreeView.cpp>
1. [DEL] Remove the code to do GPU-based POM reconstruction.
2. [MOD] Change the base class from CGlutWin to CDvrWin2 and CClipVolume. Also the methods _BeginDisplay(), _EndDisplay(), and _RenderSlab() are overloaded.
3. [ADD] Define the method _UpdateTf() to update the transfer function.
4. [ADD] In _InitFunc(), setup the transfer function and the textures and sharders for volume rendering.
5. [ADD] Define the method _TimerFunc for animation.
6. [ADD] Add a flag iIsPlottingVolume to enable/disable the volume rendering mode.
7. [ADD] Add  a flag iIsPlottingISosurface to decide whether dvr or isosurface is used.

<SOMViewer_main.cpp>
1. [ADD] Add the window for transfer function editing.
2. [ADD] Define the function _UpdateTf() to pass the latest transfer function to CSATSepDWT3DView.
3. [MOD] Change the size of the PCP view.
4. [DEL] Remove the windows 3D view, array view.
5. [DEL] Remove the sahring of display lists between the basis view and the tree view.

<SOM.h>
1. [ADD] Declate the array to store the histogram of time difference volume.
2. [ADD] Declare the array to srtore the average time difference of all grids.
3. [ADD] DEclare the matrix for fast 2nd order diff computation.

<CMakeLists.txt>
1. [ADD] Add the shader ray_integral.frag.
2. [ADD Link the libraries libdvr2, libtfw, and libclip.

<ray_integral.frag>
1. [1ST] Check in the fragment shader for volume rendering.

<SOMEmbeddingView.h>
<SOMEmbeddingView.cpp>
1. [1ST] Check in the code of the embedding view.

Revision 1.12  2012-03-15 07:49:28  leeten

[03/15/2012]
<SOM.cpp>
1. [ADD] Define the new method _ComputePathlines to reconstruct the pathlines at the given time step.
2. [MOD] When compute the time difference, divide by the number of dimension.

<SOMBasisView.cpp>
1. [MOD] Take the squared root of the accumuated 1st order.

<SOMTree.cpp>
1. [ADD] Define the method _PlacePathlines() to place the pathlines.
2. [ADD] Define the method _ClearPathlines() to clean the selected pathlines. It is also called in the desctructor().
3. [MOD] Sort the clusters in _SelectCluster().

<SOMTreeView.cpp>
1. [ADD] Define the new method _UpdateClickecSom() to update the clicked array, and use it when the mouse is clicked.
2. [ADD] Add a new GLUI control to setup the flag iIsPlottingPlacePathlines. If it is true, the placed pathlines are rendered.

<SOMViewer_main.cpp>
1. [DEL] Remove the window cSomArrayView.

<SOM.h>
1. [ADD] Declare a new method _ComputePathlines() to reconstruct the pathlines at the given time step.

<SOMTree.h>
1. [ADD] Declare a new strcutre CPathline for the pathlines, and a list to store them.
2. [ADD] Declare a new method _PlacePathlines() to place the pathlines.
3. [ADD] Declare a new method _ClearPathlines() to clean the selected pathlines.

<SOMTreeView.h>
1. [ADD] Define a new GLUI GLUI_EVENT_GENERATE_PATHLINES for the button to generate pathlines.
2. [ADD] Declare a new field fPathlineThreshold to control the threshold.
3. [ADD] Declare a new field iIsPlottingPlacePathlines to decide whether the pathlines are plotted.
4. [ADD] Declare the new method _UpdateClickecSom() to update the clicked array.

Revision 1.11  2012-03-08 21:15:14  leeten

[03/08/2012]
<SOM.cpp>
1. [ADD] Print the progress when opening the MATLAB engine.

<SOMTree.cpp>
1. [ADD] Define the method _EmbedCluster().
2. [ADD] Compute the max. of the 1st order difference norm.

<SOMTreeCluster.cpp>
1. [MOD] Change the paramete type of _Draw() from CSom* to CSomTree*.
2. [ADD] Support different coloring schemes. The cluster can be plotted by its SVD coefficients, the 1st order differecne norm, and the embedded coordinates.

<SOMTreeView.cpp>
1. [MOD] Change the paramete type of CSomTree::CCluster::_Draw() from CSom* to CSomTree*.
2. [ADD] Add GLUI radiobox to select the coloring schemes.

<SOMTree.h>
1. [ADD] Define enums as the embedding algorithms and a new field iEmbedding to index the algorithm.
2. [ADD] Define enums as the coloring schemes of clusters and a new field iColoring to index to coloring scheme.
3. [ADD] Declare the method _EmbedCluster().
4. [ADD] Declare the field fMax1stOrderDiff as the max. of the 1st order difference norm.

<SOMTreeView.h>
1. [ADD] Define a new enum ZAXIS_BEGIN_TIME_REVERSED, which means that the Z axis represente the start time in reversed order.

Revision 1.10  2012-03-04 19:34:44  leeten

[03/04/2012]
<SOMBasisView.cpp>
1. [ADD] Plot the time series of the accum diff to indicate the time span.
2. [MOD] Correct the time spans.

<SOMCluster.cpp>
1. [ADD] Compute the time series of the accum. diff in the method _Update().

<SOMTreeView.h>
<SOMTreeView.cpp>
1. [MOD] Change the type of the time spans from size_t to integer s.t. it can define the non-valid time range.

<SOMViewer_main.cpp>
1. [DEL] Remove the window SOM3DView.

<SOM.cpp>
<SOM3DView.cpp>
<SOMPcpView.cpp>
<SOMBasisView.h>
<SOMCluster.h>
<SOMTreeView.h>
1. [DEL] Remove old code segments.

Revision 1.9  2012-03-02 17:03:36  leeten

[03/02/2012]
1. [ADD] Declare the method _SetTimeRange() to specify the time range.
2. [ADD] Declare the method _SetAdvecting() to specify the advection direction.
3. [ADD] Declare 2D arrays for the SOM and the linked SOM of the clicked cluster.
4. [ADD] Define a new rendering mode to display the pathlines/streaklines in Space/Time. The scale of the time axis and the offset of the space domain can be adjusted via GLUI components.
5. [MOD] Declare a method _DisplayInSelectionMode() for the original display routine.

Revision 1.8  2012-02-15 19:12:55  leeten

[02/15/2012]
1. [DEL] Remove the fields and methods that upload the basis SOMs as textures.
2. [ADD] Declare the method _InitGL() that use OpenGL extension to setup FBOs and shaders.

Revision 1.7  2012-02-13 22:37:35  leeten

[02/13/2012]
1. [MOD] Define the structure CSomCompositeFbo in the header SOMCompositeFbo.h.
2. [ADD] Add a variable pidCompositeBasis to store the handle of the fragment shader.
3. [ADD] Define the structure CBasisTex for the basis textures.
4. [ADD] Store the f4SomRange as the ranges of the four channels.
5. [ADD] Store the additional scale as the variable fColorScale.

Revision 1.6  2012-02-09 21:40:03  leeten

[02/09/2012]
1. [ADD] Define an enum GLUI_EVENT_MAX_ID at the end so the previous one is always declared with a comma (,).
2. [ADD] Define an enum EVENT_CLICK_CLUSTER to indicate the events that a cluster has been clicked by the middle button of the mouse.
3. [ADD] Declare the method _IdleFunc() and _MouseFunc().

Revision 1.5  2012-02-01 16:55:38  leeten

[02/01/2012]
1. [DEL] Remove old deleted code segments.

Revision 1.4  2012-01-31 18:19:03  leeten

[01/31/2012]
1. [MOD] Replace the event GLUI_EVENT_ERROR by GLUI_EVENT_METRIC, GLUI_EVENT_QUADTREE_THRESHOLD, and GLUI_EVENT_CLUSTERING_THRESHOLD.
2. [ADD] Declare a new method CSomTree::CNode::BIsCovering() to check whether an given extent is totally within the given quadtree node.

Revision 1.3  2012-01-27 03:55:59  leeten

[01/26/2012]
1. [ADD] Add a new flag iIsPlottingNodes to decide whether the quadtree nodes are plotted.

Revision 1.2  2012-01-24 21:46:50  leeten

[01/24/2012]
1. [MOD] Change the parameter fErrorThreshold to fQuadtreeThreshold and fClusteringThreshold.
2. [ADD] Add flags to control the rendering of edges/clusters.

Revision 1.1  2012-01-22 23:05:21  leeten
[01/22/2012]
1. [1ST] Check in the files for the classes CSATSepDWT3DView to visualize the class CSomTree.


*/

