#include <GL/glew.h>

#include "shader.h"	
#include "libbuf3d.h"

#include "libclock.h"	

#include "SATSepDWT3DView.h"
using namespace SATSepDWT;

//////////////////// CGlutWin methods //////////////////// 
void 
C3DView::_IdleFunc()
{
	CDvrSuiteWin::_IdleFunc();	

	if( cCursor3D.bActive )
	{
		double pdTempCoord3D_obj[3];

		for(int i = 0; i < 3; i++)
			pdTempCoord3D_obj[i] = 
				cCursor3D.pdCoord_obj[i] + 
				cCursor3D.pdViewRayStep_obj[i] * (double)(iCursorY - iBeginY) / (double)piViewport[3];

		if( -1.0 <= pdTempCoord3D_obj[0] && pdTempCoord3D_obj[0] <= 1.0 &&
			-1.0 <= pdTempCoord3D_obj[1] && pdTempCoord3D_obj[1] <= 1.0 &&
			-1.0 <= pdTempCoord3D_obj[2] && pdTempCoord3D_obj[2] <= 1.0 )
		{
			for(int i = 0; i < 3; i++)
				cCursor3D.pdCoord_obj[i] = pdTempCoord3D_obj[i];

			pcQuery->i4Location = make_int4(
				(int)((double)iXDim * (cCursor3D.pdCoord_obj[0] + 1.0)/2.0),
				(int)((double)iYDim * (cCursor3D.pdCoord_obj[1] + 1.0)/2.0),
				(int)((double)iZDim * (cCursor3D.pdCoord_obj[2] + 1.0)/2.0),
				1.0);
			CGlutWin::_GlobalCB(
				IGetId(), 
				CGlutWin::CB_MANUAL, 
				EVENT_CURSOR_3D);
		}
	} 
}

void 
C3DView::_MouseFunc
	(
		int button, 
		int state, 
		int x, 
		int y
	)
{
	CDvrSuiteWin::_MouseFunc(button, state, x, y);	
	switch(button) 
	{
	case GLUT_MIDDLE_BUTTON: // pan
		cCursor3D.bActive = false;
		if( GLUT_DOWN == state ) 
		{
			switch( eModifier & ~GLUT_ACTIVE_ALT )
			{
			case 0:
				{
					TMatrix tModelViewMatrix;

					glPushMatrix();
						glLoadIdentity();
						glMultMatrixd(tModifiedModelviewMatrix);
						glGetDoublev(GL_MODELVIEW_MATRIX, tModelViewMatrix);
					glPopMatrix();

					double pdNearCoord_obj[3];
					double pdFarCoord_obj[3];

					gluUnProject(
						(double)iBeginX + 0.5, (double)iBeginY + 0.5, 0.0, 
						tModelViewMatrix, tProjectionMatrix, piViewport, 
						&pdNearCoord_obj[0], &pdNearCoord_obj[1], &pdNearCoord_obj[2]);
					
					gluUnProject(
						(double)iBeginX + 0.5, (double)iBeginY + 0.5, +1.0, 
						tModelViewMatrix, tProjectionMatrix, piViewport, 
						&pdFarCoord_obj[0], &pdFarCoord_obj[1], &pdFarCoord_obj[2]);

					double *pdViewRayOrigin_obj = pdNearCoord_obj;
					double pdViewRayDir_obj[3];
					for(int i = 0; i < 3; i++)
						pdViewRayDir_obj[i] = pdFarCoord_obj[i] - pdNearCoord_obj[i];

					double dMinT = HUGE_VAL;
					double dMaxT = -HUGE_VAL;
					bool bIntersectBox = false;
					for(int i = 0; i < 3; i++)	// direction x, y, or z
					{
						if( pdViewRayDir_obj[i] != 0.0 )
						{
							double t;
							for(int j = 0; j < 2; j++)	// +1 or -1
							{
								double w = (double)(j * 2 - 1);	// convert j to +1/-1

								t = (w - pdViewRayOrigin_obj[i]) / pdViewRayDir_obj[i];

								bool bOutsideFace = false;
								for(int k = 0; k < 3; k++)
									if( k != i )
									{
										double dCoord_obj = pdViewRayOrigin_obj[k] + t * pdViewRayDir_obj[k];
										if( -1.0 > dCoord_obj || dCoord_obj > 1.0 )
										{
											bOutsideFace = true;
											break;
										}
									}

								if( false == bOutsideFace )
								{
									bIntersectBox = true;
									dMinT = min(dMinT, t);
									dMaxT = max(dMaxT, t);
								}
							} // for j 
						} // if
					} // for i
					dMinT = max(dMinT, 0.0);
					
					if( true == bIntersectBox )
					{
						double dNrOfSteps = 
							(double)(iXDim * iXDim + iYDim * iYDim + iZDim * iZDim);
						dNrOfSteps = sqrt(dNrOfSteps);
						dNrOfSteps /= 30;	

						for(int i = 0 ; i < 3; i++ )
						{
							double dEnterCoord	= pdViewRayOrigin_obj[i] + dMinT * pdViewRayDir_obj[i];
							double dExitCoord	= pdViewRayOrigin_obj[i] + dMaxT * pdViewRayDir_obj[i];

							cCursor3D.pdViewRayStep_obj[i] = (dExitCoord - dEnterCoord) / dNrOfSteps;
							cCursor3D.pdCoord_obj[i] = dEnterCoord;
						}

						cCursor3D.bActive = true;
					}
				}
				break;
			}
		}
		break;
	} // switch(eMouseButton) 
}

void 
C3DView::
	_InitGl()
{
}

void 
C3DView::_InitFunc()
{
	CDvrSuiteWin::_InitFunc();
	_CreateTfWins();
}

void 
C3DView::
	_TimerFunc(unsigned short usEvent)
{
	CDvrSuiteWin::_TimerFunc(usEvent);	// ADD-BY-LEETEN 02/10/2013
}

void 
C3DView::_DisplayFunc()
{
	CDvrSuiteWin::_DisplayFunc();
}

void 
C3DView::_KeyboardFunc(unsigned char ubKey, int x, int y)
{
	CDvrSuiteWin::_KeyboardFunc(ubKey, x, y);	// ADD-BY-LEETEN 02/10/2013
}

void 
C3DView::_ReshapeFunc(int w, int h)
{
	CDvrSuiteWin::_ReshapeFunc(w, h);
}

void 
C3DView::_GluiFunc(unsigned short usValue)
{
	CDvrSuiteWin::_GluiFunc(usValue);	
}

//////////////////////////////////////////////////////
void 
C3DView::
	_BeginDisplay()
{
	CDvrSuiteWin::_BeginDisplay();

	glGetDoublev(GL_MODELVIEW_MATRIX, tModifiedModelviewMatrix);	
}

void 
C3DView::
	_EndDisplay()
{
	glUseProgramObjectARB(0);

	///////////////////////////////////////////////////////////
	// now plot the boxes
	glPushAttrib(
		GL_DEPTH_BUFFER_BIT |
		0);
	glDisable(GL_DEPTH_TEST);

	glPushMatrix();

	glTranslatef(-1.0f, -1.0f, -1.0f);
	glScalef(2.0f/(float)(iXDim - 1), 2.0f/(float)(iYDim - 1), 2.0f/(float)(iZDim - 1));
	pcBlockTree->_RenderBlock
	(
		CBlock::MODE_ASSIGNED 
	);
	#if	0	// DEL-BY-LEETEN 06/23/2013-BEGIN
	// ADD-BY-LEETEN 05/07/2013-BEGIN
	pcBlockTree->_RenderBlock
	(
		CBlock::MODE_SELECTED_BY_HIST
	);
	// ADD-BY-LEETEN 05/07/2013-END
	#endif	// DEL-BY-LEETEN 06/23/2013-END

	// ADD-BY-LEETEN 06/23/2013-BEGIN
	for(size_t b = 0; b < CBlock::VGetBlocksRenderedByPCP().size(); b++)
	{
		CBlock::VGetBlocksRenderedByPCP()[b]->_RenderBlock(CBlock::F4GetPCPColor());
	}
	// ADD-BY-LEETEN 06/23/2013-END

	// Plot the querying box
	if( cCursor3D.bActive )
	{
		pcQuery->i4Location.x = (int)((float)iXDim * (float)(cCursor3D.pdCoord_obj[0] + 1.0)/2.0f);
		pcQuery->i4Location.y = (int)((float)iYDim * (float)(cCursor3D.pdCoord_obj[1] + 1.0)/2.0f);
		pcQuery->i4Location.z = (int)((float)iZDim * (float)(cCursor3D.pdCoord_obj[2] + 1.0)/2.0f);
	}
	
	if( cCursor3D.bActive || pcQuery->iIsPlot )
	{
		glPushMatrix();
		glTranslatef(
			(float)pcQuery->i4Location.x, 
			(float)pcQuery->i4Location.y, 
			(float)pcQuery->i4Location.z 
			);
		glScalef(
			(float)pcQuery->i4Size.x, 
			(float)pcQuery->i4Size.y, 
			(float)pcQuery->i4Size.z); 
		glTranslatef(+0.5f, +0.5f, +0.5f);

		glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
		glutWireCube(1.0);
		glPopMatrix();
	}

	glPopMatrix();

	glPopAttrib();
		// GL_DEPTH_BUFFER_BIT |

	glPopMatrix();
}

void 
C3DView::
	_RenderSlab(
		int iSlab, int iNrOfSlabs,
		double pdModelviewMatrix[], double pdProjectionMatrix[], int piViewport[],
		double dMinX, double dMaxX, 
		double dMinY, double dMaxY, 
		double dMinZ, double dMaxZ
	)
{
	CDvrSuiteWin::_RenderSlab(
		iSlab, iNrOfSlabs,
		pdModelviewMatrix, pdProjectionMatrix, piViewport,
		dMinX, dMaxX, 
		dMinY, dMaxY, 
		dMinZ, dMaxZ);
}

//////////////////// Constructors/Destructors //////////////////// 
C3DView::
	C3DView(void)
	:
	CDvrSuiteWin()
{
	fThicknessGain = 4.0f;
	iNrOfSlices = 32;
}

C3DView::
	~C3DView(void)
{
}

