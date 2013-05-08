#pragma once

#include <string>
using namespace std;

#include <GL/glew.h>	

#include <vector_types.h>	// ADD-BY-LEETEN 02/06/2013

#include "libopt.h"
#include "liblog.h"
#include "libclock.h"	

#include "SimpleNDFile.h"
#include "SATSepDWTBlock.h"

namespace SATSepDWT 
{
	struct CBlockTree:
		public CSimpleNDFile<double>
	{
		CBlock* pcRoot;
		size_t uMaxLevel;
		vector<double> vdLevelBinMax;
		
		void
		_CompMaxProb
		(
			const int2& i2BinRange,
			void *_Reserved = NULL
		)
		{
			vdLevelBinMax.assign(this->uMaxLevel + 1, (double)-HUGE_VAL);
			vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> > vpairPDF;

			vector<size_t> vuLevel, vuLocalCoefSub, vuGlobalCoefBase, vuLocalCoefLengths;
			for(size_t l = 0; l <= this->uMaxLevel; l++)
			{
				vuLevel.assign(UGetNrOfDims(), l);
				this->_ConvertWaveletSubToLevels(vuLevel, vuGlobalCoefBase, vuLocalCoefLengths);
				size_t uNrOfLocalCoefs = 1;
				for(size_t d = 0; d < vuLocalCoefLengths.size(); d++)
					uNrOfLocalCoefs *= vuLocalCoefLengths[d];
				for(size_t lc = 0; lc < uNrOfLocalCoefs; lc++)
				{
					this->_GetCoefSparse(
						WaveletSAT::UConvertSubToIndex(vuLevel, vuDimLevels),
						lc, 
						vpairPDF);
					CBlock::_ConvertToPDF(vpairPDF);

					for(size_t bi = 0; bi < vpairPDF.size(); bi++)
					{
						WaveletSAT::typeBin uBin = vpairPDF[bi].first;
						if((int)uBin < i2BinRange.x || (int)uBin > i2BinRange.y )
							continue;
						WaveletSAT::typeWavelet dProb = vpairPDF[bi].second;
						vdLevelBinMax[l] = max(vdLevelBinMax[l], dProb);
					}
				}
			}
		}

		void
		_LoadFile
		(
			const char* szFilepath,
			void *_Reserved = NULL
		)
		{
			CSimpleNDFile<double>::_LoadFile(szFilepath);
			// Decide the max level.
			uMaxLevel = 0;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				size_t uLevel = vuDimLevels[d] - 1;
				if( !uMaxLevel )
					uMaxLevel = uLevel;
				else
					uMaxLevel = min(uMaxLevel, uLevel);
			}
			
			pcRoot = new CBlock;
			vector<size_t> vuLocalSub;
			vuLocalSub.assign(UGetNrOfDims(), 0);

			CBlock::_SetSource(this);
			const vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> > vpairCDF;
			pcRoot->_Initialize(uMaxLevel, 0, vuLocalSub, vpairCDF);

			_CompMaxProb(make_int2(0, this->UGetNrOfBins()-1));
		}
		
		virtual
		void
		_RenderHistograms(
			const CBlock::EMode eMode,
			size_t uNrOfLevelsToDisplay,
			const int2& i2BinRange,
			const float4& f4Color,
			void* _Reserved = NULL
		)
		{
			pcRoot->_RenderHistogram(
				eMode,
				uNrOfLevelsToDisplay,
				false,
				i2BinRange,
				vdLevelBinMax,
				f4Color
				);
		}

		virtual
		void
		_RenderPolylines(
			const CBlock::EMode eMode,
			const float4& f4Color,
			void* _Reserved = NULL
		)
		{
			glBegin(GL_LINES);
			pcRoot->_RenderPolyline(
				eMode,
				this->uMaxLevel,
				false,
				f4Color,
				0.0
				);
			glEnd();
		}
		virtual
		void
		_RenderBlock(
			const CBlock::EMode eMode,
			void* _Reserved = NULL
		)
		{
			pcRoot->_RenderBlock(eMode);
		}

		// ADD-BY-LEETEN 05/07/2013-BEGIN
		virtual
		void
		_FilteredPaths(
			const CBlock::EFilterAction eAction,
			const vector<float2>& vf2Filter,
			const size_t uFilterLevel,
			const float4 f4Color,
			void* _Reserved = NULL
		)
		{
			pcRoot->BFilterPaths
			(
				eAction, 
				uFilterLevel,
				vf2Filter, 
				true,
				0.0,
				f4Color
			);
		}
		// ADD-BY-LEETEN 05/07/2013-END

		virtual
		void
		_SelectByHist
		(
			const CBlock::EMode eMode,
			bool bIsAddingMode,
			const float4& f4Color,
			const int2& i2BinRange,
			const vector<float2>& vf2BinRanges,
			size_t uSelectLevel,
			void *_Reserved = NULL
		)
		{
			pcRoot->_SelectByHist(
				eMode,
				bIsAddingMode,
				f4Color,
				this->uMaxLevel,
				i2BinRange,
				vf2BinRanges,
				uSelectLevel
				);
		}

		virtual
		void
		_ResetSelection
		(
			const CBlock::EMode eResetMode,
			void *_Reserved = NULL
		)
		{
			pcRoot->_ResetSelection(eResetMode);
		}

		virtual
		~CBlockTree()
		{
			if( pcRoot )
				delete pcRoot;
		}
	};
}