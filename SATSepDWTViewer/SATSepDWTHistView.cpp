#include <GL/glew.h>

#include "shader.h"	
#include "libbuf3d.h"

#include "libclock.h"	

#include "SATSepDWTHistView.h"

// ADD-BY-LEETEN 02/03/2013-BEGIN
void
CSATSepDWTHistView::
	_Allocate(
		void *_Reserved
)
{
	CSATSepDWT::_Allocate();

	vpairCoefColors.assign(this->uNrOfCoefs, pair<bool, float4>(false, make_float4(0.0f, 0.0f, 0.0f, 1.0f)));
}

// ADD-BY-LEETEN 02/14/2013-BEGIN
void
CSATSepDWTHistView::
	_ConvertToCDF
	(
		vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >& vpairBinValue,
		void *_Reserved
	)
{
	WaveletSAT::typeWavelet dSum = (WaveletSAT::typeWavelet)0.0;
	for(size_t b = 0; b < vpairBinValue.size(); b++)
		dSum += vpairBinValue[b].second;

	// convert it to pdf first
	for(size_t b = 0; b < vpairBinValue.size(); b++)
		vpairBinValue[b].second /= dSum;

	// convert the pdf to cdf
	for(size_t b = 1; b < vpairBinValue.size(); b++)
		vpairBinValue[b].second += vpairBinValue[b-1].second;
}

WaveletSAT::typeWavelet
CSATSepDWTHistView::
	DCompCoefs
	(
		const vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >& vpairCDF1,
		const vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >& vpairCDF2,
		void *_Reserved
	)
{
	vector<WaveletSAT::typeWavelet> vCDF;
	vCDF.assign(UGetNrOfBins(), (WaveletSAT::typeWavelet)0);
	for(size_t b = 0; b < vpairCDF1.size() - 1; b++)
		for(size_t bj = vpairCDF1[b].first; bj < vpairCDF1[b+1].first; bj++)
			vCDF[bj] = vpairCDF1[b].second;
	for(size_t b = vpairCDF1[vpairCDF1.size()-1].first; b < UGetNrOfBins(); b++)
		vCDF[b] = 1.0;

	for(size_t b = 0; b < vpairCDF2.size() - 1; b++)
		for(size_t bj = vpairCDF2[b].first; bj < vpairCDF2[b+1].first; bj++)
			vCDF[bj] -= vpairCDF2[b].second;
	for(size_t b = vpairCDF2[vpairCDF2.size()-1].first; b < UGetNrOfBins(); b++)
		vCDF[b] -= 1.0;

	WaveletSAT::typeWavelet dDiff = (WaveletSAT::typeWavelet)0;
	for(size_t b = 0; b < vCDF.size(); b++)
		dDiff += vCDF[b] * vCDF[b];
	dDiff = (WaveletSAT::typeWavelet)sqrt((double)dDiff / (double)vCDF.size());

	return dDiff;
}

void
CSATSepDWTHistView::
	_CompWithChild
(
	const vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >& vpairParentCDF,
	const size_t uLevel,
	const vector<size_t>& vuLocalSub,
	void	*_Reserved
)
{
	vector<size_t> vuWaveletSub;
	vuWaveletSub.assign(UGetNrOfDims(), uLevel);
	vector<size_t> vuGlobalBase, vuLocalLengths;
	this->_ConvertWaveletSubToLevels(
		vuWaveletSub, 
		vuGlobalBase, 
		vuLocalLengths);
	vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> > vpairCDF;
	this->_GetCoefSparse(
		WaveletSAT::UConvertSubToIndex(vuWaveletSub,	vuDimLevels),
		WaveletSAT::UConvertSubToIndex(vuLocalSub,		vuLocalLengths),
		vpairCDF);
	_ConvertToCDF(vpairCDF);

	double dDiff = DCompCoefs(vpairParentCDF, vpairCDF);

	vvdDiffWithParent[uLevel]
		[WaveletSAT::UConvertSubToIndex(vuLocalSub, vuLocalLengths)] = dDiff;

	if( uLevel == uMaxLevel )
		return;

	// scan through each child, 
	size_t uNrOfChildren = (size_t)1 << UGetNrOfDims();
	vector<size_t> vuChildVol;
	vuChildVol.assign(UGetNrOfDims(), 2);
	vector<size_t> vuChildSub;
	
	for(size_t c = 0; c < uNrOfChildren; c++)
	{
		WaveletSAT::_ConvertIndexToSub(c, vuChildSub, vuChildVol);
			
		for(size_t d = 0; d < UGetNrOfDims(); d++)
			vuChildSub[d] += vuLocalSub[d] * 2;

		_CompWithChild
		(
			vpairCDF,
			uLevel + 1,
			vuChildSub
		);
	}
}

void
CSATSepDWTHistView::
	_CompMaxProb
(
	void	*_Reserved
)
{
	vector<size_t> vuLevel, vuLocalCoefSub, vuGlobalCoefBase, vuLocalCoefLengths;
	vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> > vpairCoefBinValues;
	vdLevelBinMax.assign(this->uMaxLevel, (double)-HUGE_VAL);
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
				vpairCoefBinValues);

			WaveletSAT::typeWavelet dSum = 0.0;
			for(vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >::iterator 
					ivpairCoefBinValues = vpairCoefBinValues.begin();
				ivpairCoefBinValues != vpairCoefBinValues.end();
				ivpairCoefBinValues ++)
				dSum += ivpairCoefBinValues->second;

			for(size_t bi = 0; bi < vpairCoefBinValues.size(); bi++)
			{
				WaveletSAT::typeBin uBin = vpairCoefBinValues[bi].first;
				if((int)uBin < i2BinRange.x || (int)uBin > i2BinRange.y )
					continue;
				WaveletSAT::typeWavelet dValue = vpairCoefBinValues[bi].second;
				WaveletSAT::typeWavelet dProb = dValue / dSum;
				vdLevelBinMax[l] = max(vdLevelBinMax[l], dProb);
			}
		}
	}
}
// ADD-BY-LEETEN 02/14/2013-END

// ADD-BY-LEETEN 02/06/2013-BEGIN
void
CSATSepDWTHistView::
	_LoadFile
	(
		const char* szFilepath,
		void *_Reserved
	)
{
	CSATSepDWT::_LoadFile(szFilepath);

	// decdie the max level
	// MOD-BY-LEETEN 02/03/2013-FROM:		size_t uMaxLevel = 0;
	uMaxLevel = 0;
	// MOD-BY-LEETEN 02/03/2013-END
	for(size_t d = 0; d < UGetNrOfDims(); d++)
	{
		// MOD-BY-LEETEN 02/14/2013-FROM:		size_t uLevel = vuDimLevels[d] - 2;
		size_t uLevel = vuDimLevels[d] - 1;
		// MOD-BY-LEETEN 02/14/2013-END
		if( !uMaxLevel )
			uMaxLevel = uLevel;
		else
			uMaxLevel = min(uMaxLevel, uLevel);
	}

	#if	0	// MOD-BY-LEETEN 02/14/2013-FROM:
	// now decide the range for each bin
	vvdLevelBinMax.resize(this->uMaxLevel + 1);

	vector<size_t> vuLevel, vuLocalCoefSub, vuGlobalCoefBase, vuLocalCoefLengths;
	vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> > vpairCoefBinValues;

	for(size_t l = 0; l <= this->uMaxLevel; l++)
	{
		vvdLevelBinMax[l].assign(UGetNrOfBins(), (WaveletSAT::typeWavelet)-HUGE_VAL);
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
				vpairCoefBinValues);

			WaveletSAT::typeWavelet dSum = 0.0;
			for(vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >::iterator 
					ivpairCoefBinValues = vpairCoefBinValues.begin();
				ivpairCoefBinValues != vpairCoefBinValues.end();
				ivpairCoefBinValues ++)
				dSum += ivpairCoefBinValues->second;

			for(size_t bi = 0; bi < vpairCoefBinValues.size(); bi++)
			{
				WaveletSAT::typeBin uBin = vpairCoefBinValues[bi].first;
				WaveletSAT::typeWavelet dValue = vpairCoefBinValues[bi].second;
				WaveletSAT::typeWavelet dProb = dValue / dSum;
				vvdLevelBinMax[l][uBin] = max(vvdLevelBinMax[l][uBin], dProb);
			}
		}

		// the max. of the current histogram
		WaveletSAT::typeWavelet dCurrentMax = (double)-HUGE_VAL;
		for(size_t b = UGetNrOfBins(); b > 0; b--)
		{
			dCurrentMax = max(dCurrentMax, vvdLevelBinMax[l][b - 1]);
			vvdLevelBinMax[l][b - 1] = dCurrentMax;
		}
	}	
	#else	// MOD-BY-LEETEN 02/14/2013-TO:
	i2BinRange = make_int2(0, UGetNrOfBins() - 1);
	_CompMaxProb();
	#endif	// MOD-BY-LEETEN 02/14/2013-END

	// ADD-BY-LEETEN 02/14/2013-BEGIN
	vector<size_t> vuLevel, vuLocalCoefSub, vuGlobalCoefBase, vuLocalCoefLengths;
	vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> > vpairCDF;
	this->_GetCoefSparse(0, 0, vpairCDF);
	this->_ConvertToCDF(vpairCDF);
	vector<size_t> vuLocalSub;
	vuLocalSub.assign(UGetNrOfDims(), 0);

	vvdDiffWithParent.resize(this->uMaxLevel + 1);
	for(size_t l = 0; l <= this->uMaxLevel; l++)
	{
		vector<size_t> vuWaveletSub;
		vuWaveletSub.assign(UGetNrOfDims(), l);
		this->_ConvertWaveletSubToLevels(vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths);
		size_t uNrOfLocals = 1;
		for(size_t d = 0; d < UGetNrOfDims(); d++)
			uNrOfLocals *= vuLocalCoefLengths[d];
		vvdDiffWithParent[l].resize(uNrOfLocals);
	}

	_CompWithChild(vpairCDF, 1, vuLocalSub);

	CSepDWTHeader *pcSepDWTHeader = (CSepDWTHeader*)&cSATSepDWTPCPView;
	pcSepDWTHeader->_Set(
			vuDimLengths,
			UGetNrOfBins()
			);
	cSATSepDWTPCPView._SetDiffWithParent(vvdDiffWithParent);
	cSATSepDWTPCPView.ICreate("PCP View");
	// ADD-BY-LEETEN 02/14/2013-END
}
// ADD-BY-LEETEN 02/06/2013-END

void
CSATSepDWTHistView::
	_RenderHistogram
(
	size_t uLevel,
	const vector<size_t>& vuWaveletSub,
	const vector<size_t>& vuLocalSub,
	const float4 f4Color,
	bool bIsHightLighting,
	void *_Reserved
)
{
	// fetch the coordinates and plot it 
	vector<size_t> vuGlobalCoefBase;
	vector<size_t> vuLocalCoefLengths;
	this->_ConvertWaveletSubToLevels(vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths);
	vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > > vpairBinCoefs;
	_GetCoefSparse( 
		WaveletSAT::UConvertSubToIndex(vuWaveletSub,	vuDimLevels),
		WaveletSAT::UConvertSubToIndex(vuLocalSub,		vuLocalCoefLengths),
		vpairBinCoefs );
	float fSum = 0.0f;
	for(vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > >::iterator 
			ivpairBinCoef = vpairBinCoefs.begin();
		ivpairBinCoef != vpairBinCoefs.end();
		ivpairBinCoef++)
	{
		#if	0	// DEL-BY-LEETEN 02/06/2013-BEGIN
		if( ivpairBinCoef->first < (size_t)iMinBin )
			continue;
		#endif	// DEL-BY-LEETEN 02/06/2013-END
		fSum += (float)ivpairBinCoef->second;
	}

	glPushAttrib(
		GL_LINE_BIT |
		0 );
	glLineWidth( (bIsHightLighting)?4.0:2.0 );
	glPushMatrix();
	// MOD-BY-LEETEN 02/03/2013-FROM:	glTranslatef(0.0f, (float)iMaxLevel - (float)uLevel, 0.0f);
	glTranslatef(0.0f, (float)uLevel, 0.0f);
	// MOD-BY-LEETEN 02/03/2013-END

	glColor4fv((float*)&f4Color);

	glBegin(GL_LINE_STRIP);
	int iPrevBin = 0;
	for(vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > >::iterator 
			ivpairBinCoef = vpairBinCoefs.begin();
		ivpairBinCoef != vpairBinCoefs.end();
		ivpairBinCoef++)
	{
		int iBin = (int)ivpairBinCoef->first;
		#if	0	// MOD-BY-LEETEN 02/14/2013-FROM:
		if( iBin < iMinBin )
			continue;
		#else	// MOD-BY-LEETEN 02/14/2013-TO:
		if( iBin < i2BinRange.x || iBin > i2BinRange.y )
			continue;
		#endif	// MOD-BY-LEETEN 02/14/2013-END

		if( iBin - iPrevBin > 1 )
		{
			glVertex2f((float)iPrevBin + 1, 0.0f);
			glVertex2f((float)iBin, 0.0f);
		}

		float fCount = (float)ivpairBinCoef->second;
		#if	0	// MOD-BY-LEETEN 02/06/2013-FROM:
		glVertex2f((float)iBin, fCount/fSum);
		glVertex2f((float)iBin+1, fCount/fSum);
		#else	// MOD-BY-LEETEN 02/06/2013-TO:
		float fProb = fCount/fSum;
		// MOD-BY-LEETEN 02/14/2013-FROM:		float fY = fProb / vvdLevelBinMax[uLevel][iMinBin];
		float fY = fProb / vdLevelBinMax[uLevel];
		// MOD-BY-LEETEN 02/14/2013-END
		glVertex2f((float)iBin, fY);
		glVertex2f((float)iBin+1, fY);
		#endif	// MOD-BY-LEETEN 02/06/2013-END
		iPrevBin = iBin;
	}
	if( iPrevBin != (int)UGetNrOfBins() - 1 )
	{
		glVertex2f((float)iPrevBin + 1, 0.0f);
		glVertex2f((float)UGetNrOfBins() - 1.0f, 0.0f);
	}
	glEnd();
	glPopMatrix();
	glPopAttrib();
		// GL_LINE_BIT |
}
// ADD-BY-LEETEN 02/03/2013-END

void
CSATSepDWTHistView::_RenderBlock
(
	size_t uLevel,
	const vector<size_t>& vuWaveletSub,
	const vector<size_t>& vuLocalSub,
	// MOD-BY-LEETEN 02/03/2013-FROM:	const float pfColor[]
	const float4 f4Color,
	bool	bIsHighLighting,
	bool	bIsNotRecursive,	// ADD-BY-LEETEN 02/11/2013
	void	*_Reserved
	// MOD-BY-LEETEN 02/03/2013-END
)
{
	if( (int)uLevel > iMaxLevel )
		return;

	#if	0	// MOD-BY-LEETEN 02/03/2013-FROM:
	// fetch the coordinates and plot it 
	vector<size_t> vuGlobalCoefBase;
	vector<size_t> vuLocalCoefLengths;
	this->_ConvertWaveletSubToLevels(vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths);
	vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > > vpairBinCoefs;
	_GetCoefSparse( 
		WaveletSAT::UConvertSubToIndex(vuWaveletSub,	vuDimLevels),
		WaveletSAT::UConvertSubToIndex(vuLocalSub,		vuLocalCoefLengths),
		vpairBinCoefs );
	float fSum = 0.0f;
	for(vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > >::iterator 
			ivpairBinCoef = vpairBinCoefs.begin();
		ivpairBinCoef != vpairBinCoefs.end();
		ivpairBinCoef++)
	{
		if( ivpairBinCoef->first < (size_t)iMinBin )
			continue;
		fSum += (float)ivpairBinCoef->second;
	}

	glPushMatrix();
	glTranslatef(0.0f, (float)iMaxLevel - (float)uLevel, 0.0f);
	glColor4fv(pfColor);
	glBegin(GL_LINE_STRIP);
	int iPrevBin = 0;
	for(vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > >::iterator 
			ivpairBinCoef = vpairBinCoefs.begin();
		ivpairBinCoef != vpairBinCoefs.end();
		ivpairBinCoef++)
	{
		int iBin = (int)ivpairBinCoef->first;
		if( iBin < iMinBin )
			continue;

		if( iBin - iPrevBin > 1 )
		{
			glVertex2f((float)iPrevBin + 1, 0.0f);
			glVertex2f((float)iBin, 0.0f);
		}

		float fCount = (float)ivpairBinCoef->second;
		glVertex2f((float)iBin, fCount/fSum);
		glVertex2f((float)iBin+1, fCount/fSum);
		iPrevBin = iBin;
	}
	if( iPrevBin != (int)UGetNrOfBins() - 1 )
	{
		glVertex2f((float)iPrevBin + 1, 0.0f);
		glVertex2f((float)UGetNrOfBins() - 1.0f, 0.0f);
	}
	glEnd();
	glPopMatrix();
	#else	// MOD-BY-LEETEN 02/03/2013-TO:
	float4 f4NextColor;

	_RenderHistogram
	(
		uLevel,
		vuWaveletSub,
		vuLocalSub,
		f4Color,
		bIsHighLighting
	);
	#endif	// MOD-BY-LEETEN 02/03/2013-END
	// ADD-BY-LEETEN 02/11/2013-BEGIN
	if(bIsNotRecursive)
		return;
	// ADD-BY-LEETEN 02/11/2013-END
	// ADD-BY-LEETEN 02/14/2013-BEGIN
	if( (int)uLevel >= iMaxLevel )
		return;
	// ADD-BY-LEETEN 02/14/2013-END
	///////////////////////////////////////////
	// now plot the children

	// first, find the global base in the next level
	vector<size_t> vuNextWaveletSub;
	vuNextWaveletSub.resize(UGetNrOfDims());
	for(size_t d = 0; d < UGetNrOfDims(); d++)
		vuNextWaveletSub[d] = vuWaveletSub[d] + 1;
		
	// scan through each child, 
	size_t uNrOfChildren = (!uLevel)?1:((size_t)1 << UGetNrOfDims());
	vector<size_t> vuChildVol;
	vuChildVol.assign(UGetNrOfDims(), 2);
	vector<size_t> vuChildSub;
	vuChildSub.resize(UGetNrOfDims());
	
	// ADD-BY-LEETEN 02/03/2013-BEGIN
	///////////////////////////////////////////////////////////////////////
	// if this children has its own color, save it to a queue first
	// display other children that do not have its own color
	// display the children that has its own color

	#if	0	// DEL-BY-LEETEN 02/10/2013-BEGIN
	// queue of the children of its own color
	vector<size_t> vuChildrenWithColor;
	#endif	// DEL-BY-LEETEN 02/10/2013-END
	vector<size_t> vuNextGlobalCoefBase, vuNextLocalCoefLengths;
	this->_ConvertWaveletSubToLevels(vuNextWaveletSub, vuNextGlobalCoefBase, vuNextLocalCoefLengths);
	// ADD-BY-LEETEN 02/03/2013-END

	for(size_t c = 0; c < uNrOfChildren; c++)
	{
		WaveletSAT::_ConvertIndexToSub(c, vuChildSub, vuChildVol);
			
		for(size_t d = 0; d < UGetNrOfDims(); d++)
			vuChildSub[d] += vuLocalSub[d] * 2;

		// ADD-BY-LEETEN 02/03/2013-BEGIN
		vector<size_t> vuGlobalSub;
		vuGlobalSub.resize(UGetNrOfDims());
		for(size_t d = 0; d < vuGlobalSub.size(); d++)
			vuGlobalSub[d] = vuNextGlobalCoefBase[d] + vuChildSub[d];
		size_t uGlobal = WaveletSAT::UConvertSubToIndex(vuGlobalSub, vuCoefLengths);

		if( this->vpairCoefColors[uGlobal].first )
		{
			vuChildrenWithColor.push_back(uGlobal);
			continue;
		}
		// ADD-BY-LEETEN 02/03/2013-END

		// otherwise, use the passed colors
		_RenderBlock(
			uLevel + 1, 
			vuNextWaveletSub, 
			vuChildSub, 
			// MOD-BY-LEETEN 02/03/2013-FROM:			pfColor);
			f4Color,
			bIsHighLighting);
			// MOD-BY-LEETEN 02/03/2013-END
	}
	#if	0	// DEL-BY-LEETEN 02/10/2013-BEGIN
	// ADD-BY-LEETEN 02/03/2013-BEGIN
	vector<size_t> 
		// MOD-BY-LEETEN 02/06/2013-FROM:		vuBlobalCoefBase,	// not used
		vuGlobalCoefBase,
		// MOD-BY-LEETEN 02/06/2013-END
		vuLocalCoefLengths;	// not used
	for(size_t c = 0; c < vuChildrenWithColor.size(); c++)
	{
		this->_ConvertIndexToLevels(
			vuChildrenWithColor[c], 
			vuNextWaveletSub, 
			vuChildSub, 
			vuGlobalCoefBase,		// not used
			vuLocalCoefLengths);	// not used
			
		// if this child has its own cluster, use the color of this cluster

		// otherwise, use the passed colors
		_RenderBlock(
			uLevel + 1, 
			vuNextWaveletSub, 
			vuChildSub, 
			this->vpairCoefColors[vuChildrenWithColor[c]].second,
			false);

		// ADD-BY-LEETEN 02/06/2013-BEGIN
		// append the blocks with its own color to vpairBlockColors
		float4 f4Left, f4Size;
		float* pfLeft = (float*)&f4Left.x;
		float* pfSize = (float*)&f4Size.x;
		for(size_t d = 0; d < UGetNrOfDims(); d++)
		{
			float fWaveletLength = (float)vuCoefLengths[d] / (float)vuLocalCoefLengths[d];
			pfLeft[d] = (float)vuChildSub[d] * fWaveletLength;
			pfSize[d] = (float)fWaveletLength;
		}
		for(size_t d = UGetNrOfDims(); d < 4; d++)
		{
			pfLeft[d] = 0.0f;
			pfSize[d] = 1.0f;
		}
		vpairBlockColors.push_back(
				make_pair<pair<float4, float4>, float4>
				(
					make_pair<float4, float4>(f4Left, f4Size), 
					this->vpairCoefColors[vuChildrenWithColor[c]].second
				)
			);
		// ADD-BY-LEETEN 02/06/2013-END
	}
	// ADD-BY-LEETEN 02/03/2013-END
	#endif	// DEL-BY-LEETEN 02/10/2013-END
}

//////////////////// CGlutWin methods //////////////////// 
void 
CSATSepDWTHistView::_IdleFunc()
{
	// ADD-BY-LEETEN 02/06/2013-BEGIN
	if( iIsPlottingBoxs )
	{
		CGlutWin::_GlobalCB(
			IGetId(), 
			CGlutWin::CB_MANUAL, 
			EVENT_PLOTTING_BOX, 
			iIsPlottingBoxs,
			(vector< pairBlockColor >*)&vpairBlockColors,
			NULL);		
	}
	// ADD-BY-LEETEN 02/06/2013-END
}

void 
CSATSepDWTHistView::
	_MotionFunc
	(
		int x, 
		int y
	)
{
#if	0	
	if( GLUT_LEFT_BUTTON == iButton )
	{
		if( pcEditingBox )
		{
			pcEditingBox->ppairCorners[1].first  = ...;
			pcEditingBox->ppairCorners[1].second = ... ;

		}
	}
#endif
}

void 
CSATSepDWTHistView::_MouseFunc
	(
		int button, 
		int state, 
		int x, 
		int y
	)
{
}

void 
CSATSepDWTHistView::
	_InitGl()
{
}

void 
CSATSepDWTHistView::_InitFunc()
{
	_KeepUpdateOn();
	
	/////////////////////////////////////////////
	// set up GLUI
	GLUI *pcGlui = PCGetGluiSubwin();
	GLUI_Spinner* pcSpinner_MaxLevel = pcGlui->add_spinner("Max Level", GLUI_SPINNER_INT, &iMaxLevel);
	#if	0	// DEL-BY-LEETEN 02/06/2013-BEGIN
		// decdie the max level
		// MOD-BY-LEETEN 02/03/2013-FROM:		size_t uMaxLevel = 0;
		uMaxLevel = 0;
		// MOD-BY-LEETEN 02/03/2013-END
		for(size_t d = 0; d < UGetNrOfDims(); d++)
		{
			size_t uLevel = vuDimLevels[d] - 2;
			if( !uMaxLevel )
				uMaxLevel = uLevel;
			else
				uMaxLevel = min(uMaxLevel, uLevel);
		}
	#endif	// DEL-BY-LEETEN 02/06/2013-END	
		pcSpinner_MaxLevel->set_int_limits(0, (int)uMaxLevel);
	#if	0	// MOD-BY-LEETEN 02/14/2013-FROM:
	GLUI_Spinner* pcSpinner_MinBin = pcGlui->add_spinner("Min Bin", GLUI_SPINNER_INT, &iMinBin);
		pcSpinner_MinBin->set_int_limits(0, UGetNrOfBins());
	#else	// MOD-BY-LEETEN 02/14/2013-TO:
	GLUI_Spinner* pcSpinner_MaxBin = pcGlui->add_spinner(
			"Max Bin", GLUI_SPINNER_INT, &i2BinRange.y, 
			IAddWid(GLUI_EVENT_BIN_RANGE), CGlutWin::_GluiCB_static);
		pcSpinner_MaxBin->set_int_limits(0, UGetNrOfBins() - 1);
	GLUI_Spinner* pcSpinner_MinBin = pcGlui->add_spinner(
			"Min Bin", GLUI_SPINNER_INT, &i2BinRange.x,
			IAddWid(GLUI_EVENT_BIN_RANGE), CGlutWin::_GluiCB_static);
		pcSpinner_MinBin->set_int_limits(0,	 UGetNrOfBins() - 1);
	#endif	// MOD-BY-LEETEN 02/14/2013-END
	// ADD-BY-LEETEN 02/11/2013-BEGIN
	{
		GLUI_Panel *pcPanel_Color = pcGlui->add_panel("Default Color");
		static char* pszChannels[] = {"R", "G", "B", "A"};
		float *pfColor = &f4DefaultColor.x;
		for(int c = 0; c < sizeof(pszChannels)/sizeof(pszChannels[0]); c++)
		{
			GLUI_Spinner* pcSpinner = pcGlui->add_spinner_to_panel(pcPanel_Color, pszChannels[c], GLUI_SPINNER_FLOAT, &pfColor[c]);
			pcSpinner->set_float_limits(0.0f, 1.0f);
		}
	}
	// ADD-BY-LEETEN 02/11/2013-END

	// MOD-BY-LEETEN 02/11/2013-FROM:	pcGlui->add_checkbox("Plotting Box?", &iIsPlottingBoxs);	// ADD-BY-LEETEN 02/06/2013
	pcGlui->add_checkbox("Plot Boxes?", &iIsPlottingBoxs,
				IAddWid(GLUI_EVENT_PLOT_BOXES), CGlutWin::_GluiCB_static);
	pcGlui->add_checkbox("Show Max. Prob.?", &iIsShowingMaxProb);
	pcGlui->add_checkbox("Not Recursive?", &iIsNotRecursive);
	// MOD-BY-LEETEN 02/11/2013-END

	// ADD-BY-LEETEN 02/03/2013-BEGIN
	/////////////////////////////////////////////////////////
	// setup the color editor
	cColorEditor.AddGlui(
		(CGlutWin*)this, 
		pcGlui, 
		NULL,
		uMaxLevel);
	// ADD-BY-LEETEN 02/03/2013-END

	/////////////////////////////////////////////////////////
	// setup the cluster editor
	#if	0	// MOD-BY-LEETEN 02/03/2013-FROM:
	GLUI_Panel *pcPanel_Cluster = pcGlui->add_panel("Cluster");
	{
		GLUI_Panel *pcPanel_Color = pcGlui->add_panel_to_panel(pcPanel_Cluster, "Color");
		static char* pszChannels[] = {"R", "G", "B", "A"};
		for(int c = 0; c < sizeof(pszChannels)/sizeof(pszChannels[0]); c++)
		{
			GLUI_Spinner* pcSpinner = pcGlui->add_spinner_to_panel(pcPanel_Color , pszChannels[c], GLUI_SPINNER_FLOAT, &cEditing.cCluster.pfColor[c]);
			pcSpinner->set_float_limits(0.0f, 1.0f);
		}
		GLUI_Spinner* pcSpinner_Level = pcGlui->add_spinner_to_panel(pcPanel_Cluster, "Level",	
				GLUI_SPINNER_INT, &cEditing.iLevel, 
				IAddWid(GLUI_EVENT_CLUSTER_EDITING), CGlutWin::_GluiCB_static);
			pcSpinner_Level->set_int_limits(0, uMaxLevel);
		GLUI_Spinner* pcSpinner_ID = pcGlui->add_spinner_to_panel(pcPanel_Cluster, "ID",		
				GLUI_SPINNER_INT, &cEditing.iID,	
				IAddWid(GLUI_EVENT_CLUSTER_EDITING), CGlutWin::_GluiCB_static);
			pcSpinner_ID->set_int_limits(0, CCluster::NR_OF_CLUSTERS_PER_LEVEL - 1);
		GLUI_Spinner* pcSpinner_Bin = pcGlui->add_spinner_to_panel(pcPanel_Cluster, "Bin",		
				GLUI_SPINNER_INT, &cEditing.iBin,
				IAddWid(GLUI_EVENT_CLUSTER_EDITING), CGlutWin::_GluiCB_static);
			pcSpinner_Bin->set_int_limits(0, UGetNrOfBins() - 1);

		GLUI_Spinner* pcSpinner_Prob;
		pcSpinner_Prob = pcGlui->add_spinner_to_panel(pcPanel_Cluster, "Upper",		
				GLUI_SPINNER_FLOAT, &cEditing.f2Prob.y,
				IAddWid(GLUI_EVENT_CLUSTER_EDITING), CGlutWin::_GluiCB_static);
			pcSpinner_Prob->set_float_limits(0.0f, 1.0f);
		pcSpinner_Prob = pcGlui->add_spinner_to_panel(pcPanel_Cluster, "Lower",		
				GLUI_SPINNER_FLOAT, &cEditing.f2Prob.x,
				IAddWid(GLUI_EVENT_CLUSTER_EDITING), CGlutWin::_GluiCB_static);
			pcSpinner_Prob->set_float_limits(0.0f, 1.0f);

		pcGlui->add_checkbox_to_panel(pcPanel_Cluster, "Editing?", &iIsEditingCluster);
	}

	cEditing.cCluster = CCluster(UGetNrOfBins());

	// allocate the clusters
	vvcClusters.resize(uMaxLevel + 1);
	for(size_t l = 0; l < vvcClusters.size(); l++)
	{
		vvcClusters[l].resize(CCluster::NR_OF_CLUSTERS_PER_LEVEL);
		for(size_t c = 0; c < vvcClusters[l].size(); c++)
			vvcClusters[l][c] = CCluster(UGetNrOfBins());
	}
	#else	// MOD-BY-LEETEN 02/03/2013-TO:
	cClusterEditor.AddGlui(
		(CGlutWin*)this, 
		pcGlui, 
		NULL,
		uMaxLevel,
		UGetNrOfBins());
	#endif	// MOD-BY-LEETEN 02/03/2013-END
}

void 
CSATSepDWTHistView::
	_TimerFunc(unsigned short usEvent)
{
}

void 
CSATSepDWTHistView::_DisplayFunc()
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// convert the screen space to [0 ... #Bins] x [0 ... #Levels]
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef(-1.0, -1.0, 0.0f);
	glScalef(2.0, 2.0, 1.0f);
	#if	0	// MOD-BY-LEETEN 02/14/2013-FROM:
	glScalef(1.0f/((float)UGetNrOfBins() - 1.0f - (float)iMinBin), 1.0f/(float)(iMaxLevel + 1), 1.0f);
	glTranslatef(-(float)iMinBin, 0.0f, 0.0f);
	#else	// MOD-BY-LEETEN 02/14/2013-TO:
	glScalef(1.0f/(float)(i2BinRange.y - i2BinRange.x), 1.0f/(float)(iMaxLevel + 1), 1.0f);
	glTranslatef(-(float)i2BinRange.x, 0.0f, 0.0f);
	#endif	// MOD-BY-LEETEN 02/14/2013-END

	////////////////////////////////////////////////////////////////
	// plot the axis
	for(size_t l = 0; l <= (size_t)iMaxLevel; l++)
	{
		glPushMatrix();
		// MOD-BY-LEETEN 02/11/2013-FROM:		glTranslatef(0.0f, (float)iMaxLevel - l, 0.0f);
		glTranslatef(0.0f, (float)l, 0.0f);
		// MOD-BY-LEETEN 02/11/2013-END

		// plot the axis
		glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
		glBegin(GL_LINES);
		glVertex2i(0, 0);
		glVertex2i((int)UGetNrOfBins(), 0);
		glEnd();

		// ADD-BY-LEETEN 02/11/2013-BEGIN
		if( iIsShowingMaxProb )
		{
			char szMaxProb[1024];
			#if	0	// MOD-BY-LEETEN 02/14/2013-FROM:
			sprintf(szMaxProb, "Prob = %f", vvdLevelBinMax[l][iMinBin]);
			_DrawString3D(szMaxProb, (float)iMinBin, 0.88f);		
			#else	// MOD-BY-LEETEN 02/14/2013-TO:
			sprintf(szMaxProb, "Prob = %f", vdLevelBinMax[l]);
			_DrawString3D(szMaxProb, (float)i2BinRange.x, 0.88f);		
			#endif	// MOD-BY-LEETEN 02/14/2013-END
		}
		// ADD-BY-LEETEN 02/11/2013-END
		glPopMatrix();
	}

	////////////////////////////////////////////////////////////////// 
	// plot block histograms
	// ADD-BY-LEETEN 02/11/2013-BEGIN
	glPushAttrib(
		GL_COLOR_BUFFER_BIT |
		0 );
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	// ADD-BY-LEETEN 02/11/2013-END
	// MOD-BY-LEETEN 02/03/2013-FROM:	static float pfDefaultColor[4] = {0.0f, 0.0f, 0.0f, 1.0f};
	// DEL-BY-LEETEN 02/11/2013:	float4 f4DefaultColor = make_float4(0.0f, 0.0f, 0.0f, 1.0f);
	// MOD-BY-LEETEN 02/03/2013-END
	vector<size_t> vuWaveletSub;
	vuWaveletSub.assign(UGetNrOfDims(), 0);
	vector<size_t> vuLocalSub;
	vuLocalSub.assign(UGetNrOfDims(), 0);
	// ADD-BY-LEETEN 02/06/2013-BEGIN
	// reset the block colors
	vpairBlockColors.clear();
	// ADD-BY-LEETEN 02/06/2013-END
	vuChildrenWithColor.clear();	// ADD-BY-LEETEN 02/10/2013
	_RenderBlock
		(
			0, 
			vuWaveletSub,
			vuLocalSub,
			// MOD-BY-LEETEN 02/03/2013-FROM:	pfDefaultColor
			f4DefaultColor,
			false
			// MOD-BY-LEETEN 02/03/2013-END
		);

	// ADD-BY-LEETEN 02/10/2013-BEGIN
	// render all blocks with its own colors from all levels here
	vector<size_t> 
		vuLevel, 
		vuLocalCoefSubs,
		vuGlobalCoefBase,
		vuLocalCoefLengths;	// not used
	for(size_t c = 0; c < vuChildrenWithColor.size(); c++)
	{
		size_t uGlobal = vuChildrenWithColor[c];
		this->_ConvertIndexToLevels(
			uGlobal, 
			vuLevel, 
			vuLocalCoefSubs, 
			vuGlobalCoefBase,		// not used
			vuLocalCoefLengths);	// not used
			
		// if this child has its own cluster, use the color of this cluster

		// otherwise, use the passed colors
		_RenderBlock(
			vuLevel[0], 
			vuLevel, 
			vuLocalCoefSubs, 
			this->vpairCoefColors[uGlobal].second,
			// MOD-BY-LEETEN 02/11/2013-FROM:			false);
			false,	// do not highlight
			(iIsNotRecursive)?true:false	// not recursive
			);
			// MOD-BY-LEETEN 02/11/2013-END

		// append the blocks with its own color to vpairBlockColors
		float4 f4Left, f4Size;
		float* pfLeft = (float*)&f4Left.x;
		float* pfSize = (float*)&f4Size.x;
		for(size_t d = 0; d < UGetNrOfDims(); d++)
		{
			float fWaveletLength = (float)vuCoefLengths[d] / (float)vuLocalCoefLengths[d];
			pfLeft[d] = (float)vuLocalCoefSubs[d] * fWaveletLength;
			pfSize[d] = (float)fWaveletLength;
		}
		for(size_t d = UGetNrOfDims(); d < 4; d++)
		{
			pfLeft[d] = 0.0f;
			pfSize[d] = 1.0f;
		}
		vpairBlockColors.push_back(
				make_pair<pair<float4, float4>, float4>
				(
					make_pair<float4, float4>(f4Left, f4Size), 
					this->vpairCoefColors[vuChildrenWithColor[c]].second
				)
			);
	}
	// ADD-BY-LEETEN 02/10/2013-END
	// ADD-BY-LEETEN 02/11/2013-BEGIN
	glPopAttrib();
		// GL_COLOR_BUFFER_BIT |
	// ADD-BY-LEETEN 02/11/2013-END
	// ADD-BY-LEETEN 02/03/2013-BEGIN
	if( cColorEditor.iIsActive )
	{
		// ADD-BY-LEETEN 02/11/2013-BEGIN
		// plot the current histogram as a dashed lin
		glPushAttrib(
			GL_LINE_BIT | 
			0);
		glLineWidth(4.0f);
		glLineStipple(4, 0xCCCC);
		glEnable(GL_LINE_STIPPLE);
		// ADD-BY-LEETEN 02/11/2013-END
		bool bIsValid;
		vector<size_t> vuWaveletSub, vuLocalSub, vuGlobalSub;
		size_t uGlobal;
		_GetColorEditor(bIsValid, vuWaveletSub, vuLocalSub, vuGlobalSub, uGlobal);
		if(bIsValid)
		{
			_RenderBlock
			(
				cColorEditor.iLevel, 
				vuWaveletSub,
				vuLocalSub,
				cColorEditor.f4Color,
				true
			);

			// ADD-BY-LEETEN 02/06/2013-BEGIN
			vector<size_t> vuGlobalCoefBase, vuLocalCoefLengths;
			this->_ConvertWaveletSubToLevels(
				vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths);

			// append the blocks with its own color to vpairBlockColors
			float4 f4Left, f4Size;
			float* pfLeft = (float*)&f4Left.x;
			float* pfSize = (float*)&f4Size.x;
			for(size_t d = 0; d < UGetNrOfDims(); d++)
			{
				float fWaveletLength = (float)vuCoefLengths[d] / (float)vuLocalCoefLengths[d];
				pfLeft[d] = (float)vuLocalSub[d] * fWaveletLength;
				pfSize[d] = (float)fWaveletLength;
			}
			for(size_t d = UGetNrOfDims(); d < 4; d++)
			{
				pfLeft[d] = 0.0f;
				pfSize[d] = 1.0f;
			}
			vpairBlockColors.push_back(
					make_pair<pair<float4, float4>, float4>
					(
						make_pair<float4, float4>(f4Left, f4Size), 
						cColorEditor.f4Color
					)
				);
			// ADD-BY-LEETEN 02/06/2013-END
		}
		// ADD-BY-LEETEN 02/11/2013-BEGIN
		glPopAttrib();
			// GL_LINE_BIT | 
		// ADD-BY-LEETEN 02/11/2013-END
	}
	// ADD-BY-LEETEN 02/03/2013-END

	////////////////////////////////////////////////////////////////
	// plot the editing cluster
	// MOD-BY-LEETEN 02/03/2013-FROM:	if( iIsEditingCluster )
	if( cClusterEditor.iIsActive )
	// MOD-BY-LEETEN 02/03/2013-END
	{
		glPushAttrib(
			GL_LINE_BIT | 
			0);
		glLineWidth(4.0f);
		// MOD-BY-LEETEN 02/11/2013-FROM:		glLineStipple(4, 0x0F0F);
		glLineStipple(4, 0xCCCC);
		// MOD-BY-LEETEN 02/11/2013-END
		glEnable(GL_LINE_STIPPLE);
		glPushMatrix();
		#if	0	// MOD-BY-LEETEN 02/03/2013-FROM:
		glTranslatef(0.0f, (float)iMaxLevel - this->cEditing.iLevel, 0.0f);
		// plot the lines
		glColor4f(this->cEditing.cCluster.pfColor[0], this->cEditing.cCluster.pfColor[1], this->cEditing.cCluster.pfColor[2], 0.1f);
		#else	// MOD-BY-LEETEN 02/03/2013-TO:
		glTranslatef(0.0f, (float)cClusterEditor.iLevel, 0.0f);
		// plot the lines
		glColor4f(cClusterEditor.f4Color.x, cClusterEditor.f4Color.y, cClusterEditor.f4Color.z, 0.1f);
		#endif	// MOD-BY-LEETEN 02/03/2013-END
		for(size_t i = 0; i < 2; i++)
		{
			glBegin(GL_LINE_STRIP);
			// MOD-BY-LEETEN 02/14/2013-FROM:			for(size_t b = this->iMinBin; b < this->UGetNrOfBins(); b++)
			for(size_t b = i2BinRange.x; b <= i2BinRange.y; b++)
			// MOD-BY-LEETEN 02/14/2013-END
			{
				// MOD-BY-LEETEN 02/03/2013-FROM:	float fProb = (i)?this->cEditing.cCluster.vf2BinRanges[b].x:this->cEditing.cCluster.vf2BinRanges[b].y;
				float fProb = (i)?this->cClusterEditor.vf2BinRanges[b].x:cClusterEditor.vf2BinRanges[b].y;
				// ADD-BY-LEETEN 02/06/2013-BEGIN
				#if	0	// MOD-BY-LEETEN 02/14/2013-FROM:
				fProb = max(min(fProb, (float)vvdLevelBinMax[cClusterEditor.iLevel][this->iMinBin]), 0.0f);
				fProb /= vvdLevelBinMax[cClusterEditor.iLevel][this->iMinBin];
				#else	// MOD-BY-LEETEN 02/14/2013-TO:
				fProb = max(min(fProb, (float)vdLevelBinMax[cClusterEditor.iLevel]), 0.0f);
				fProb /= vdLevelBinMax[cClusterEditor.iLevel];
				#endif	// MOD-BY-LEETEN 02/14/2013-END
				// ADD-BY-LEETEN 02/06/2013-END
				// MOD-BY-LEETEN 02/03/2013-END
				glVertex2f((float)b, fProb);
				glVertex2f((float)b + 1.0f, fProb);
			}
			glEnd();
		}
		// ADD-BY-LEETEN 02/03/2013-BEGIN
		// plot a rectangle to highlight the current bin
		glDisable(GL_LINE_STIPPLE);
		float l = (float)cClusterEditor.iBin;
		float r = l + 1.0f;
		float b = cClusterEditor.f2Prob.x;
		float t = cClusterEditor.f2Prob.y;
		// ADD-BY-LEETEN 02/06/2013-BEGIN
		#if	0	// MOD-BY-LEETEN 02/14/2013-FROM:
		b = max(min(b, (float)vvdLevelBinMax[cClusterEditor.iLevel][this->iMinBin]), 0.0f) 
			/ vvdLevelBinMax[cClusterEditor.iLevel][this->iMinBin];
		t = max(min(t, (float)vvdLevelBinMax[cClusterEditor.iLevel][this->iMinBin]), 0.0f)
			/ vvdLevelBinMax[cClusterEditor.iLevel][this->iMinBin];
		#else	// MOD-BY-LEETEN 02/14/2013-TO:
		b = max(min(b, (float)vdLevelBinMax[cClusterEditor.iLevel]), 0.0f) / vdLevelBinMax[cClusterEditor.iLevel];
		t = max(min(t, (float)vdLevelBinMax[cClusterEditor.iLevel]), 0.0f) / vdLevelBinMax[cClusterEditor.iLevel];
		#endif	// MOD-BY-LEETEN 02/14/2013-END
		// ADD-BY-LEETEN 02/06/2013-END
		glBegin(GL_LINE_STRIP);
		glVertex2f(l, b);
		glVertex2f(r, b);
		glVertex2f(r, t);
		glVertex2f(l, t);
		glVertex2f(l, b);
		glEnd();
		// ADD-BY-LEETEN 02/03/2013-END
		glPopMatrix();
		glPopAttrib(
			// GL_LINE_BIT
		);
	}

	// restore the coordinate system back
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void 
CSATSepDWTHistView::_KeyboardFunc(unsigned char ubKey, int x, int y)
{
}

void 
CSATSepDWTHistView::_ReshapeFunc(int w, int h)
{
}

// ADD-BY-LEETEN 02/03/2013-BEGIN
void
CSATSepDWTHistView::
	_GetColorEditor(
		bool& bIsValid,
		vector<size_t>& vuWaveletSub,
		vector<size_t>& vuLocalSub,
		vector<size_t>& vuGlobalSub,
		size_t& uGlobal,
		void* _Reserved
	)
{
	bIsValid = true;
	if( !cColorEditor.iIsActive )
		bIsValid = false;

	if( cColorEditor.iLevel > (int)this->uMaxLevel || cColorEditor.iLevel < 0 )
		bIsValid = false;

	if( !bIsValid )
		return;

	vuWaveletSub.assign(UGetNrOfDims(), cColorEditor.iLevel);

	vector<size_t> vuGlobalCoefBase, vuLocalCoefLengths;
	this->_ConvertWaveletSubToLevels(vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths);

	size_t uNrOfBlocks = 1;
	for(size_t d = 0; d < vuLocalCoefLengths.size(); d++)
		uNrOfBlocks *= vuLocalCoefLengths[d];

	if( cColorEditor.iBlock >= (int)uNrOfBlocks || cColorEditor.iBlock < 0 )
		bIsValid = false;

	if( !bIsValid )
		return;

	// now convert the local ID to local subscript
	WaveletSAT::_ConvertIndexToSub(cColorEditor.iBlock, vuLocalSub, vuLocalCoefLengths);

	vuGlobalSub.resize(UGetNrOfDims());
	for(size_t d = 0; d < vuLocalCoefLengths.size(); d++)
		vuGlobalSub[d] = vuGlobalCoefBase[d] + vuLocalSub[d];

	uGlobal = WaveletSAT::UConvertSubToIndex(vuGlobalSub, vuCoefLengths);
}
// ADD-BY-LEETEN 02/03/2013-END

void 
CSATSepDWTHistView::_GluiFunc(unsigned short usValue)
{
	switch(usValue)
	{
	// ADD-BY-LEETEN 02/14/2013-BEGIN
	case GLUI_EVENT_BIN_RANGE:
	{
		_CompMaxProb();
	}	break;
	// ADD-BY-LEETEN 02/14/2013-END
	// ADD-BY-LEETEN 02/11/2013-BEGIN
	case GLUI_EVENT_PLOT_BOXES:
	{
		if( !iIsPlottingBoxs )
		{
			CGlutWin::_GlobalCB(
				IGetId(), 
				CGlutWin::CB_MANUAL, 
				EVENT_PLOTTING_BOX, 
				iIsPlottingBoxs,
				(vector< pairBlockColor >*)&vpairBlockColors,
				NULL);		
		}
	} break;
	// ADD-BY-LEETEN 02/11/2013-END
	// ADD-BY-LEETEN 02/03/2013-BEGIN
	case GLUI_EVENT_COLOR_ASSIGN:
	{
		bool bIsValid;
		vector<size_t> vuWaveletSub, vuLocalSub, vuGlobalSub;
		size_t uGlobal;
		_GetColorEditor(bIsValid, vuWaveletSub, vuLocalSub, vuGlobalSub, uGlobal);
		if( !bIsValid )
			return;

		vpairCoefColors[uGlobal].first = true;
		vpairCoefColors[uGlobal].second = cColorEditor.f4Color;
	} break;

	case GLUI_EVENT_COLOR_RESET:
	{
		bool bIsValid;
		vector<size_t> vuWaveletSub, vuLocalSub, vuGlobalSub;
		size_t uGlobal;
		_GetColorEditor(bIsValid, vuWaveletSub, vuLocalSub, vuGlobalSub, uGlobal);
		if( !bIsValid )
			return;
		vpairCoefColors[uGlobal].first = false;
	} break;
	// ADD-BY-LEETEN 02/03/2013-END

	// ADD-BY-LEETEN 02/03/2013-BEGIN
	case GLUI_EVENT_CLUSTER_RESET_PROB:
		cClusterEditor.f2Prob = make_float2(0.0f, 1.0f);
		cClusterEditor.vf2BinRanges[cClusterEditor.iBin] = cClusterEditor.f2Prob;	// ADD-BY-LEETEN 02/11/2013
		break;

	case GLUI_EVENT_CLUSTER_ASSIGN:
	case GLUI_EVENT_CLUSTER_RESET:
	{
		vector<size_t> vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths;
		vuWaveletSub.assign(UGetNrOfDims(), cClusterEditor.iLevel);
		size_t uWavelet = WaveletSAT::UConvertSubToIndex(vuWaveletSub,	vuDimLevels);
		this->_ConvertWaveletSubToLevels(vuWaveletSub, vuGlobalCoefBase, vuLocalCoefLengths);
		size_t uNrOfChildren = 1;
		for(size_t d = 0; d < vuLocalCoefLengths.size(); d++)
			uNrOfChildren *= vuLocalCoefLengths[d];
		vector<size_t> vuLocalSub;
		for(size_t c = 0; c < uNrOfChildren; c++)
		{
			WaveletSAT::_ConvertIndexToSub(c, vuLocalSub, vuLocalCoefLengths);
			vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > > vpairBinCoefs;
			_GetCoefSparse( 
				uWavelet,
				c, 
				vpairBinCoefs);
			float fSum = 0.0f;
			for(size_t b = 0; b < vpairBinCoefs.size(); b++)
				fSum += vpairBinCoefs[b].second;

			vector< pair< WaveletSAT::typeBin, WaveletSAT::typeWavelet > >::iterator ivpairBinCoef = vpairBinCoefs.begin();
			bool bIsIn = true;
			// MOD-BY-LEETEN 02/14/2013-FROM:			for(size_t b = this->iMinBin; b < this->UGetNrOfBins(); b++)
			for(size_t b = i2BinRange.x; b <= i2BinRange.y; b++)
			// MOD-BY-LEETEN 02/14/2013-END
			{
				// update the bin
				float fProb = 0.0f;
				if( vpairBinCoefs.end() == ivpairBinCoef )
					fProb = 0.0f;
				else
				{
					if( b < ivpairBinCoef->first )
						fProb = 0.0f;
					else
					{
						fProb = ivpairBinCoef->second / fSum;
						ivpairBinCoef++;
					}
				}

				if( fProb < cClusterEditor.vf2BinRanges[b].x || fProb > cClusterEditor.vf2BinRanges[b].y )
				{
					bIsIn = false;
					break;
				}
			}
			if( !bIsIn )
				continue;

			vector<size_t> vuGlobalSub;
			vuGlobalSub.resize(UGetNrOfDims());
			for(size_t d = 0; d < vuGlobalSub.size(); d++)
				vuGlobalSub[d]  = vuGlobalCoefBase[d] + vuLocalSub[d];
			size_t uGlobal = WaveletSAT::UConvertSubToIndex(vuGlobalSub, vuCoefLengths);
			switch(usValue)
			{
				case GLUI_EVENT_CLUSTER_ASSIGN:
					this->vpairCoefColors[uGlobal].first = true;
					this->vpairCoefColors[uGlobal].second = cClusterEditor.f4Color;
					break;
				case GLUI_EVENT_CLUSTER_RESET:
					this->vpairCoefColors[uGlobal].first = false;
					break;
			}
		}
	} break;
	// ADD-BY-LEETEN 02/03/2013-END

	case GLUI_EVENT_CLUSTER_EDITING:
	{
		// MOD-BY-LEETEN 02/03/2013-FROM:		cEditing.cCluster.vf2BinRanges[cEditing.iBin] = cEditing.f2Prob;
		cClusterEditor.vf2BinRanges[cClusterEditor.iBin] = cClusterEditor.f2Prob;
		// MOD-BY-LEETEN 02/03/2013-END
	} break;
	}
}

//////////////////// Constructors/Destructors //////////////////// 
CSATSepDWTHistView::
	CSATSepDWTHistView(void)
{
	f4DefaultColor = make_float4(0.0f, 0.0f, 0.0f, 0.1f);	// ADD-BY-LEETEN 02/11/2013

	// add a panel for the UI control
	_AddGluiSubwin(GLUI_SUBWINDOW_LEFT);
}

CSATSepDWTHistView::
	~CSATSepDWTHistView(void)
{
}

