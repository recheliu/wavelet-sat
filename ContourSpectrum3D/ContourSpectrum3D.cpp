#include <vector>
using namespace std;

#include <math.h>
#include <assert.h>
#include <stdlib.h> 
#include "libclock.h"
#include "libopt.h"

#include "NrrdIO.h"
#include "SimpleND.h"
#include "contourspectrum.h"

Nrrd *nin;
#if	!WITH_DOUBLE_COEF	
	typedef float typeData;
#else	// #if	!WITH_DOUBLE_COEF
	typedef double typeData;
#endif	// #if	!WITH_DOUBLE_COEF
typeData	dValueMin = (typeData)HUGE_VAL;
typeData	dValueMax = (typeData)-HUGE_VAL;
vector<typeData> vdData;

//! Convert the volume to an array of double type
template<typename T>
void 
_ConvertVolume
(
	const Nrrd *nin,
	typeData& dValueMin,
	typeData& dValueMax,
	vector<typeData>& vdData
)
{	
	T *data = (T*)nin->data;	

	// search for the range
	dValueMin = (typeData)HUGE_VAL;
	dValueMax = (typeData)-HUGE_VAL;
	for(int v = 0,	z = 0; z < (int)nin->axis[2].size; z++)
		for(int		y = 0; y < (int)nin->axis[1].size; y++)
			for(int x = 0; x < (int)nin->axis[0].size; x++, v++)
			{
				typeData dValue = (typeData)data[v];
				vdData.push_back(dValue);
				dValueMin = min(dValueMin, dValue);
				dValueMax = max(dValueMax, dValue);
			}
	LOG_VAR(dValueMin);
	LOG_VAR(dValueMax);
}

void
_ReadVolume
(
	const char* szPathFilename
)
{
	/* create a nrrd; at this point this is just an empty container */
	nin = nrrdNew();

	/* tell nrrdLoad to only read the header, not the data */
	NrrdIoState *nio = nrrdIoStateNew();
	nrrdIoStateSet(nio, nrrdIoStateSkipData, AIR_TRUE);

	/* read in the nrrd from file */
	if (nrrdLoad(nin, szPathFilename, nio)) {
		char *err = biffGetDone(NRRD);
		LOG_ERROR(fprintf(stderr, "%s", err));
		free(err);
		return;
	}

	/* we're done with the nrrdIoState, this sets it to NULL */
	nio = nrrdIoStateNix(nio);

	LOG_VAR(nrrdElementNumber(nin));
	LOG_VAR(nrrdElementSize(nin));
	nin->data = calloc(nrrdElementNumber(nin), nrrdElementSize(nin));

	if (nrrdLoad(nin, szPathFilename, NULL)) {
		char *err = biffGetDone(NRRD);
		LOG_ERROR(fprintf(stderr, "%s", err));
		free(err);
		return;
	}

	switch(nin->type)
	{
	case nrrdTypeChar:	_ConvertVolume<char>(nin, dValueMin, dValueMax, vdData);		break;
	case nrrdTypeUChar:	_ConvertVolume<unsigned char>(nin, dValueMin, dValueMax, vdData);	break;
	case nrrdTypeShort:	_ConvertVolume<short>(nin, dValueMin, dValueMax, vdData);		break;
	case nrrdTypeUShort:	_ConvertVolume<unsigned short>(nin, dValueMin, dValueMax, vdData);	break;
	case nrrdTypeInt:	_ConvertVolume<int>(nin, dValueMin, dValueMax, vdData);			break;
	case nrrdTypeUInt:	_ConvertVolume<unsigned int>(nin, dValueMin, dValueMax, vdData);	break;
	case nrrdTypeFloat:	_ConvertVolume<float>(nin, dValueMin, dValueMax, vdData);		break;

	default:
		break;
	}
}

int
main(int argn, char* argv[])
{
	_OPTInit();			// initialize the option parser

	char *szVolFilePath = NULL;
	_OPTAddStringVector(
		"--vol-filepath", 1,
		&szVolFilePath, szVolFilePath);

	int iNrOfBins = 8;	// iValueMax;
	_OPTAddIntegerVector(
		"--n-bins", 1,
		&iNrOfBins, iNrOfBins);

	bool bIsOptParsed = BOPTParse(argv, argn, 1);

	assert(bIsOptParsed);
	assert(szVolFilePath);

	LOG_VAR(szVolFilePath);
	_ReadVolume(szVolFilePath);

	size_t uNrOfDims = (size_t)nin->dim;
	size_t uMaxLevel = 0;
	size_t uWinSize = 1;
	size_t uNrOfBins = (size_t)iNrOfBins;

	vector<size_t> vuDimLengths;
	for(size_t d = 0; d < uNrOfDims; d++)
	{
		size_t uDimLength = (size_t)nin->axis[d].size;
		vuDimLengths.push_back( uDimLength );
	}

	vector<pair<double, glm::dvec4> > vCorners;
	vCorners.resize(8);

	vector<double> vdContourSpectrum;
	vdContourSpectrum.assign(uNrOfBins, 0.0);
	vector<double> vdCellHistogram;
	vdCellHistogram.assign(uNrOfBins, 0.0);
	size_t uFirstCellBin = uNrOfBins;
	vector<double> vdBinEdges;
	vdBinEdges.resize(uNrOfBins+1);
	typeData dBinInterval = (dValueMax - dValueMin)/(double)uNrOfBins;
	for(size_t b = 0; b < vdBinEdges.size(); b++)
		vdBinEdges[b] = (double)b * dBinInterval;

	for(size_t			z = 0; z < nin->axis[2].size - 1; z++)
		for(size_t		y = 0; y < nin->axis[1].size - 1; y++)
			for(size_t	x = 0; x < nin->axis[0].size - 1; x++)
			{
				size_t uTetraOrientation = 0;
				if( 1 == x % 2 )
					uTetraOrientation = 1 - uTetraOrientation;
				if( 1 == y % 2 )
					uTetraOrientation = 1 - uTetraOrientation;
				if( 1 == z % 2 )
					uTetraOrientation = 1 - uTetraOrientation;
				for(size_t j = 0,	zo = 0; zo < 2; zo++)
					for(size_t		yo = 0; yo < 2; yo++)
						for(size_t	xo = 0; xo < 2; xo++, j++)
						{
							vector<size_t> vuSubs;
							vuSubs.resize(3);
							vuSubs[0] = x + xo;
							vuSubs[1] = y + yo;
							vuSubs[2] = z + zo;
							vCorners[j].first = vdData[WaveletSAT::UConvertSubToIndex(vuSubs, vuDimLengths)];	
							vCorners[j].second = glm::dvec4((double)x + xo, (double)y + yo, (double)z + zo, 1.0);
						}

				ContourSpectrum::_ComputeFor3DCell
				(
					vdBinEdges,
					vCorners,
					uFirstCellBin,
					vdCellHistogram,
					uTetraOrientation
				);
				for(size_t b = uFirstCellBin; b < vdCellHistogram.size(); b++)
				{
					if( !vdCellHistogram[b] )
						break;
					vdContourSpectrum[b] += vdCellHistogram[b];
				}
			}

	vector<double> vdHistogram;
	vdHistogram.assign(uNrOfBins, 0.0);
	for(size_t	v = 0,	z = 0; z < nin->axis[2].size; z++)
		for(size_t		y = 0; y < nin->axis[1].size; y++)
			for(size_t	x = 0; x < nin->axis[0].size; x++, v++)
			{
				size_t uBin = (size_t)floor((vdData[v] - dValueMin)/dBinInterval);
				uBin = min(uBin, uNrOfBins - 1);
				vdHistogram[uBin] += 1.0;
			}

	for(size_t b = 0; b < vdContourSpectrum.size(); b++)
		printf("B: %f, %f\n", vdHistogram[b], vdContourSpectrum[b]);

	return 0;
}

/***********************************************************
Copyright (c) 2013 Teng-Yok Lee

See the file LICENSE.txt for copying permission.
***********************************************************/
