#pragma once

#include <string>
using namespace std;

#include <GL/glew.h>	
#include <GL/glut.h>

#include <vector_types.h>	
#include <vector_functions.h>	

#include "libopt.h"
#include "liblog.h"
#include "libclock.h"	

#include "SimpleNDFile.h"

namespace SATSepDWT 
{
	struct CBlock 
	{
		enum EMode {
			MODE_NONE,
			MODE_ASSIGNED = 0x0001,
			MODE_SELECTED_BY_PCP = 0x0002,
			MODE_SELECTED_BY_HIST = 0x0004,
		};

		typedef double typeData;
		static CSimpleNDFile<typeData>* pcSimpleND;

		static 
		void
		_SetSource
		( 
			CSimpleNDFile<typeData>* _pcSimpleND,
			void *_Reserved = NULL)
		{
			pcSimpleND = _pcSimpleND;
		}

		static 
		CSimpleNDFile<typeData>& 
		CGetSource
		(
			void *_Reserved = NULL)
		{
			return *pcSimpleND;
		}

		static 
		WaveletSAT::typeWavelet 
		DCompCoefs
		(
			const vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >& vpairCDF1,
			const vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >& vpairCDF2,
			void *_Reserved = NULL
		)
		{
			if( vpairCDF1.empty() || vpairCDF2.empty() )
				return ( vpairCDF1.empty() && vpairCDF2.empty() )?0.0:1.0;
			vector<WaveletSAT::typeWavelet> vCDF;
			vCDF.assign(CGetSource().UGetNrOfBins(), 0.0);
			for(size_t b = 0; b < vpairCDF1.size() - 1; b++)
				for(size_t bj = vpairCDF1[b].first; bj < vpairCDF1[b+1].first; bj++)
					vCDF[bj] = vpairCDF1[b].second;
			for(size_t b = vpairCDF1[vpairCDF1.size()-1].first; b < CGetSource().UGetNrOfBins(); b++)
				vCDF[b] = 1.0;

			for(size_t b = 0; b < vpairCDF2.size() - 1; b++)
				for(size_t bj = vpairCDF2[b].first; bj < vpairCDF2[b+1].first; bj++)
					vCDF[bj] -= vpairCDF2[b].second;
			for(size_t b = vpairCDF2[vpairCDF2.size()-1].first; b < CGetSource().UGetNrOfBins(); b++)
				vCDF[b] -= 1.0;

			WaveletSAT::typeWavelet dDiff = 0.0;
			for(size_t b = 0; b < vCDF.size(); b++)
				dDiff += vCDF[b] * vCDF[b];
			dDiff = sqrt(dDiff / (WaveletSAT::typeWavelet)vCDF.size());

			return dDiff;
		}

		static
		void
		_ConvertToPDF
		(
			vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >& vpairBinValue,
			void *_Reserved = NULL
		)
		{
			WaveletSAT::typeWavelet dSum = 0.0;
			for(size_t b = 0; b < vpairBinValue.size(); b++)
				dSum += vpairBinValue[b].second;

			// convert it to pdf first
			for(size_t b = 0; b < vpairBinValue.size(); b++)
				vpairBinValue[b].second /= dSum;
		}

		static
		void
		_ConvertToCDF
		(
			vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >& vpairBinValue,
			void *_Reserved = NULL
		)
		{
			_ConvertToPDF(vpairBinValue);

			// convert the pdf to cdf
			for(size_t b = 1; b < vpairBinValue.size(); b++)
				vpairBinValue[b].second += vpairBinValue[b-1].second;
		}
		/////////////////////////////////////////////////////////////////////////

		unsigned int eMode;
		float4 f4Color;
		size_t uLevel;

		pair<float4, float4> pairExtent;

		// Use these two values to retrieve the coefficients
		size_t uWavelet;	
		size_t uLocalCoef;

		vector<CBlock *> vpcChildren;

		WaveletSAT::typeWavelet dDistFromParent;
		
		bool
		BIsModeMatch( 
			const EMode eRef,
			void *_Reserved = NULL
			)
		{
			if( MODE_NONE == eRef )
			{
				return (eMode == eRef);
			}
			else
			{
				return (eMode & eRef > 0);
			}
		}

		void
		_ModifyMode(
			const EMode eNewMode, 
			const bool bIsAddingMode,
			const float4& f4NewColor,
			void *_Reserved = NULL
		)
		{
			if( MODE_NONE == eNewMode )
			{
				if( bIsAddingMode )
					eMode = MODE_NONE;
			}
			else
			{
				if( MODE_ASSIGNED == eNewMode && bIsAddingMode )
				{
					this->f4Color = f4NewColor;
				}

				if( bIsAddingMode )
				{
					eMode |= eNewMode;
				}
				else
					eMode &= ~eNewMode;
			}
		}

		void
		_Initialize(
			const size_t uMaxLevel,
			const size_t uLevel, 
			const vector<size_t>& vuSub,
			const vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >&  vpairParentCDF,
			void *_Reserved = NULL)
		{
			vector<size_t> vuWaveletSub;
			vuWaveletSub.assign(CGetSource().UGetNrOfDims(), uLevel);
			uWavelet = WaveletSAT::UConvertSubToIndex(vuWaveletSub, CGetSource().VGetDimLevels());

			// Compute the wavelet ID.
			vector<size_t> vuGlobalCoefBase, vuLocalCoefLengths;
			CGetSource()._ConvertWaveletSubToLevels(vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths);

			// Compute the local coef. ID.
			this->uLocalCoef = WaveletSAT::UConvertSubToIndex(vuSub, vuLocalCoefLengths);

			// Fetch the CDF.
			vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> > vpairCDF;
			CGetSource()._GetCoefSparse(uWavelet, uLocalCoef, vpairCDF);
			_ConvertToCDF(vpairCDF);
			this->dDistFromParent = ( 0 == uLevel )?0.0:this->dDistFromParent = DCompCoefs(vpairParentCDF, vpairCDF);

			// Decide the extent
			const vector<size_t>& vuCoefLengths = CGetSource().VGetCoefLengths();
			float* pfFirst = &this->pairExtent.first.x;
			float* pfSecond = &this->pairExtent.second.x;
			for(size_t d = 0; d < CGetSource().UGetNrOfDims(); d++)
			{
				pfFirst[d] = (float)vuCoefLengths[d] * (float)vuSub[d] / (float)vuLocalCoefLengths[d];
				pfSecond[d] = (float)vuCoefLengths[d] * (float)(vuSub[d] + 1) / (float)vuLocalCoefLengths[d];
			}

			this->uLevel= uLevel;

			if( uLevel == uMaxLevel )
				return;

			// create the children
			size_t uNrOfChildren;
			vector<size_t> vuChildrenDimLengths;
			if( !uLevel )
			{
				uNrOfChildren = 1;
				vuChildrenDimLengths.assign(CGetSource().UGetNrOfDims(), 1);
			}
			else
			{
				uNrOfChildren = (size_t)1 << CGetSource().UGetNrOfDims();
				vuChildrenDimLengths.assign(CGetSource().UGetNrOfDims(), 2);
			}

			vpcChildren.assign(uNrOfChildren, NULL);
			vector<size_t> vuLocalSub;
			for(size_t c = 0; c < uNrOfChildren; c++) 
			{
				WaveletSAT::_ConvertIndexToSub(c, vuLocalSub, vuChildrenDimLengths);
				for(size_t d = 0; d < CGetSource().UGetNrOfDims(); d++)
					vuLocalSub[d] += vuSub[d] * 2;
				vpcChildren[c] = new CBlock();
				vpcChildren[c]->_Initialize(
					uMaxLevel, 
					uLevel + 1, 
					vuLocalSub, 
					vpairCDF);
			}
		}		

		void
		_SelectByPcp
		(
			void *_Reserved = NULL
		)
		{

		}

		void
		_ResetSelection
		(
			const EMode eResetMode,
			void *_Reserved = NULL
		)
		{
			_ModifyMode(eResetMode, false, f4Color);
			for(size_t c = 0; c < vpcChildren.size(); c++)
				vpcChildren[c]->_ResetSelection(eResetMode);
		}

		void
		_SelectByHist
		(
			const EMode eNewMode,
			const bool bIsAddingMode,
			const float4& f4Color,
			const size_t uMaxLevel,
			const int2& i2BinRange,
			const vector<float2>& vf2BinRanges,
			size_t uSelectLevel,
			void *_Reserved = NULL
		)
		{
			if( uLevel == uSelectLevel )
			{
				vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> > vpairPDF;
				CGetSource()._GetCoefSparse(uWavelet, uLocalCoef, vpairPDF);
				_ConvertToPDF(vpairPDF);
				bool bIsIn = true;
				for(size_t pdfi = 0, b = i2BinRange.x; b <= i2BinRange.y; b++)
				{
					float fProb = 0.0f;
					for(;vpairPDF[pdfi].first < b && pdfi<vpairPDF.size();pdfi++)
						;
					if(vpairPDF[pdfi].first == b)
						fProb = (float)vpairPDF[pdfi].second;

					if( fProb < vf2BinRanges[b].x || fProb > vf2BinRanges[b].y )
					{
						bIsIn = false;
						break;
					}
				}
				if( bIsIn ) 
				{
					_ModifyMode(eNewMode, bIsAddingMode, f4Color);
				}
			}
			if( uLevel < uSelectLevel || uLevel < uMaxLevel )
				for(size_t c = 0; c < vpcChildren.size(); c++)
					vpcChildren[c]->_SelectByHist
					(
						eNewMode,
						bIsAddingMode,
						f4Color,
						uMaxLevel,
						i2BinRange,
						vf2BinRanges,
						uSelectLevel
					);
		}

		void
		_RenderPolyline
		(
			const EMode eRenderMode,
			const size_t uNrOfLevelsToDisplay,
			const bool bIsAnscentMatched,
			const float4& f4DefaultColor,
			const double dParentDist,
			void *_Reserved = NULL
		)
		{
			float4 f4NewColor = f4DefaultColor;
			bool bIsRender = false;
			bool bIsMatch = bIsAnscentMatched;
			bool bIsRecursive = true;
			_SetupRendering
			(
				eRenderMode,
				f4DefaultColor,
				bIsAnscentMatched,
				uNrOfLevelsToDisplay,

				f4NewColor,
				bIsRender,
				bIsMatch,
				bIsRecursive
			);

			if( bIsRender )
			{
				glColor4fv(&f4NewColor.x);
				glVertex2d((double)uLevel - 1,	dParentDist);
				glVertex2d((double)uLevel,		this->dDistFromParent);
			} 

			if( bIsRecursive )
				for(size_t c = 0; c < this->vpcChildren.size(); c++) 
					vpcChildren[c]->_RenderPolyline
					(
						eRenderMode,
						uNrOfLevelsToDisplay, 
						bIsMatch,
						f4NewColor,
						this->dDistFromParent
					);
		}

		void
		_RenderBlock
		(
			const EMode eRenderMode,
			void *_Reserved = NULL
		)
		{
			float4 f4NewColor = f4Color;
			bool bIsRender = false;
			if( MODE_ASSIGNED == eRenderMode && this->BIsModeMatch(eRenderMode) )
				bIsRender = true;

			if( bIsRender )
			{
				glColor4fv(&f4NewColor.x);
				glPushMatrix();
				const float4& f4Left = this->pairExtent.first;
				glTranslatef(f4Left.x, f4Left.y, f4Left.z);
				const float4& f4Right = this->pairExtent.second;
				glScalef(f4Right.x - f4Left.x, f4Right.y - f4Left.y, f4Right.z - f4Left.z);
				glTranslatef(+0.5f, +0.5f, +0.5f);
				glColor4fv(&f4NewColor.x);
				glutWireCube(1.0);
				glPopMatrix();
			} 

			for(size_t c = 0; c < this->vpcChildren.size(); c++) 
				vpcChildren[c]->_RenderBlock
				(
					eRenderMode
				);
		}

		void
		_SetupRendering
		(
			const EMode eRenderMode,
			const float4& f4DefaultColor,
			const bool& bIsAnscentMatched,
			const size_t uMaxLevel,

			float4& f4NewColor,
			bool& bIsRender,
			bool& bIsMatch,
			bool& bIsRecursive,

			void *_Reserved = NULL
		)
		{
			f4NewColor = f4DefaultColor;
			bIsRender = false;
			bIsMatch = bIsAnscentMatched;
			bIsRecursive = true;
			if( MODE_NONE == eRenderMode )
			{
				if( BIsModeMatch(eRenderMode) )
					bIsRender = true;
			}
			else
			{
				if( bIsAnscentMatched )
					bIsRender = true;
				else
				if( BIsModeMatch(eRenderMode) )
				{
					bIsRender = true;
					bIsMatch = true;
					if( MODE_ASSIGNED == eRenderMode )
						f4NewColor = f4Color;
				}
			}
			if( uLevel >= uMaxLevel )
				bIsRecursive = false;
			else
			{
				if(MODE_NONE == eRenderMode && !BIsModeMatch(eRenderMode) )
					bIsRecursive = false;
			}
		}

		void
		_RenderHistogram
		(
			const EMode eRenderMode,
			const size_t uNrOfLevelsToDisplay,
			const bool bIsAnscentMatched,
			const int2& i2BinRange,
			const vector<double>& vdLevelBinMax,
			const float4& f4DefaultColor,
			void *_Reserved = NULL
		)
		{
			float4 f4NewColor = f4DefaultColor;
			bool bIsRender = false;
			bool bIsMatch = bIsAnscentMatched;
			bool bIsRecursive = true;
			_SetupRendering
			(
				eRenderMode,
				f4DefaultColor,
				bIsAnscentMatched,
				uNrOfLevelsToDisplay,

				f4NewColor,
				bIsRender,
				bIsMatch,
				bIsRecursive
			);

			if( bIsRender )
			{
				vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > > vpairPDF;
				CGetSource()._GetCoefSparse(this->uWavelet, this->uLocalCoef, vpairPDF );
				_ConvertToPDF(vpairPDF);

				glPushAttrib(
					GL_LINE_BIT |
					0 );
				// glLineWidth( (bIsHightLighting)?4.0:2.0 );
				glPushMatrix();
				glTranslatef(0.0f, (float)uLevel, 0.0f);

				glColor4fv((float*)&f4NewColor);
				glBegin(GL_LINE_STRIP);
				int iPrevBin = 0;
				for(vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > >::iterator 
						ivpairPDF = vpairPDF.begin();
					ivpairPDF != vpairPDF.end();
					ivpairPDF++)
				{
					int iBin = (int)ivpairPDF->first;
					if( iBin < i2BinRange.x || iBin > i2BinRange.y )
						continue;

					if( iBin - iPrevBin > 1 )
					{
						glVertex2f((float)iPrevBin + 1, 0.0f);
						glVertex2f((float)iBin, 0.0f);
					}

					float fProb = (float)ivpairPDF->second;
					float fY = fProb / vdLevelBinMax[uLevel];
					// fY *= (float)dWeight;

					glVertex2f((float)iBin, fY);
					glVertex2f((float)iBin+1, fY);
					iPrevBin = iBin;
				}
				if( iPrevBin != (int)CGetSource().UGetNrOfBins() - 1 )
				{
					glVertex2f((float)iPrevBin + 1, 0.0f);
					glVertex2f((float)CGetSource().UGetNrOfBins() - 1.0f, 0.0f);
				}
				glEnd();
				glPopMatrix();
				glPopAttrib();
					// GL_LINE_BIT |
			}

			if( bIsRecursive )
				for(size_t c = 0; c < vpcChildren.size(); c++)
				{
					vpcChildren[c]->_RenderHistogram(
						eRenderMode,
						uNrOfLevelsToDisplay,
						bIsMatch,
						i2BinRange,
						vdLevelBinMax,
						f4NewColor
					);
				}
		}

		CBlock()
		{
			eMode = MODE_NONE;
		}

		virtual
		~CBlock()
		{
			for(size_t c = 0; c < this->vpcChildren.size(); c++)
				delete this->vpcChildren[c];
		}
	};
}