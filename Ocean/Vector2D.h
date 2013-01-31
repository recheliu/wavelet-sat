#pragma once

#include "SimpleND.h"

template<
	typename DT = WaveletSAT::typeBin,		//!< Type of the data. Here the data is the bin
	typename ST = WaveletSAT::typeSum,		//!< Type of the sum
	typename BT = WaveletSAT::typeBin,		//!< Type of the bin
	typename WT = WaveletSAT::typeWavelet	//!< Type of the wavelet coefficientsd
>
class CVector2D:
		virtual public CSimpleND<DT, ST, BT, WT>
{
	typedef pair<DT, DT> typeVec;

	enum {
		BIN_MISSING_VALUE,
		BIN_ZERO_LENGTH,
		NR_OF_DEFAULT_BINS
	};

	enum{
		NR_OF_COMPS = 2,
	};

	vector<float> pvfComps[NR_OF_COMPS];
	float pfMissingValues[NR_OF_COMPS];
	vector<DT> vuBins;

	size_t uDepth;

	size_t uNrOfDataBins;
public:
	enum EParameter {
		PARAMETER_BEGIN = 0x1000,
		NR_OF_BINS,
		DEPTH,			
		PARAMETER_END
	};

	virtual	
	void
	_SetInteger(
		int eName,
		long lValue,
		void* _Reserved = NULL
	)
	{
		switch(eName)
		{
		case NR_OF_BINS:
			// here one more bin is allocated for the missing value
			uNrOfDataBins = (size_t)lValue;	
			// MOD-BY-LEETEN 01/31/2013-FROM:			uNrOfBins = uNrOfDataBins + NR_OF_DEFAULT_BINS;	
			this->uNrOfBins = uNrOfDataBins + NR_OF_DEFAULT_BINS;	
			// MOD-BY-LEETEN 01/31/2013-END
			break;

		case DEPTH:
			uDepth = (size_t)lValue;
			break;
		}
		CSimpleND<DT, ST, BT, WT>::_SetInteger(eName, lValue);
	}

	virtual
	void
	_ReadVectorComponent
	(
		size_t uComp,
		const char *szVecDirPath,
		const char *szVecFileName,
		const char *szVecVarName,
		void* _Reserved = NULL
	)
	{
		if( uComp >= NR_OF_COMPS )
		{
			LOG_ERROR(cerr<<"Invalid uComp");
			return;
		}

		char szVecFilePath[NC_MAX_NAME];
		sprintf(szVecFilePath, "%s/%s", szVecDirPath, szVecFileName);

		int ncId;
		ASSERT_NETCDF(nc_open(
    		szVecFilePath,
    		NC_NOWRITE,
    		&ncId));

		int varComp;
		ASSERT_NETCDF(nc_inq_varid(
			ncId, 
			szVecVarName, 
			&varComp) );

		ASSERT_NETCDF(nc_get_att_float(
			ncId,
			varComp,
			"missing_value",
			&pfMissingValues[uComp]));

		int iNrOfDims;
		ASSERT_NETCDF(nc_inq_varndims(
			ncId, 
			varComp,
			&iNrOfDims));

		int pdimIDs[NC_MAX_DIMS];
		ASSERT_NETCDF(nc_inq_vardimid(
			ncId, 
			varComp,
			pdimIDs));

		size_t uDataSize = 1;
		vector<size_t> vuDimLengths;
		vuDimLengths.clear();
		for(size_t d = 0; d < (size_t) NR_OF_COMPS; d++)
		{
			size_t uDimLen;
			ASSERT_NETCDF(nc_inq_dimlen(
				ncId, 
				pdimIDs[iNrOfDims - 1 - d], 
				&uDimLen));
			uDataSize *= uDimLen;

			vuDimLengths.push_back(uDimLen);
		}

		// set up the dim length for this slice
		size_t puStarts[NC_MAX_DIMS];
		size_t puCounts[NC_MAX_DIMS];
		for(size_t d = 0; d < (size_t)(NR_OF_COMPS + 1); d++)
			if( 0 ==  d)
			{
				puStarts[d] = uDepth;
				puCounts[d] = 1;
			}
			else
			{
				puStarts[d] = 0;
				puCounts[d] = vuDimLengths[NR_OF_COMPS - d];
			}

		// load the data
		vector<float>& vfComp = pvfComps[uComp];
		vfComp.resize(uDataSize);
		ASSERT_NETCDF(nc_get_vara(
			ncId, 
			varComp, 
			&puStarts[0],
			&puCounts[0], 
			vfComp.data()));

		ASSERT_NETCDF(nc_close(ncId));

		// if the componenet ID is larger than 1, this method can directly return
		if( uComp )
			return;

		// now set up the encoder
		// MOD-BY-LEETEN 01/31/2013-FROM:		_Set(vuDimLengths, uNrOfBins);
		_Set(vuDimLengths, this->uNrOfBins);
		// MOD-BY-LEETEN 01/31/2013-END
	}

	virtual
	size_t 
	UGetBin(
		double dU, 
		double dV,
		void* _Reserved = NULL
		)
	{
		size_t uBin = 0;
		if( dU == pfMissingValues[0] && dV == pfMissingValues[1] )
		{
			uBin = BIN_MISSING_VALUE;
		}
		else
		{
			double dRadiusSqaured = dU * dU + dV * dV;
			double dAngle_radian = atan2(dV, dU);
			if( 0.0 == dRadiusSqaured )
			{
				uBin = BIN_ZERO_LENGTH;
			}
			else
			{
				uBin = (size_t)floor((double)uNrOfDataBins * (dAngle_radian + M_PI) / (2.0 * M_PI));
				uBin = min(uBin, uNrOfDataBins - 1);
				uBin += NR_OF_DEFAULT_BINS; 
			}
		}
		return uBin;
	}

	// load the vector field and then setup the encoder
	virtual
	void
	_Load
	(
		const char *szVecDirPath,
		const char *szUVecFileName,
		const char *szUVecVarName,
		const char *szVVecFileName,
		const char *szVVecVarName,
		void* _Reserved = NULL
	)
	{
		_ReadVectorComponent(0, szVecDirPath, szUVecFileName, szUVecVarName);
		_ReadVectorComponent(1, szVecDirPath, szVVecFileName, szVVecVarName);
		const size_t& uDataSize = this->uDataSize; // ADD-BY-LEETEN 01/31/2013
		vuBins.assign(uDataSize, 0);
		for(size_t d = 0; d < uDataSize; d++)
			vuBins[d] = (DT)UGetBin(pvfComps[0][d], pvfComps[1][d]);

		// MOD-BY-LEETEN 01/31/2013-FROM:		_Allocate();
		this->_Allocate();
		// MOD-BY-LEETEN 01/31/2013-END
	}

	virtual 
	void 
	_MapValueToBins
	(
		const vector<size_t>& vuPos,
		const DT& value, 
		vector< pair<BT, ST> >& vpBins,
		void *_Reserved = NULL
	)
	{
		BT uBin = (BT)value;
		vpBins.clear();
		vpBins.push_back(pair<BT, ST>((BT)uBin, (ST)1));
	}

	////////////////////////////////////////////////////////////
	virtual
	void
	_Encode
	(
		void *_Reserved = NULL
	)
	{
	  // ADD-BY-LEETEN 01/31/2013-BEGIN
	  const size_t& uDataSize = this->uDataSize;
	  const vector<size_t>& vuDimLengths = vuDimLengths;
	  // ADD-BY-LEETEN 01/31/2013-END

		for(size_t d = 0; d < uDataSize; d++)
		{
			vector<size_t> vuPos;
			WaveletSAT::_ConvertIndexToSub(d, vuPos, vuDimLengths);
			this->_Update(vuPos, vuBins[d]);
		}
	  // MOD-BY-LEETEN 01/31/2013-FROM: 		_Finalize();
		this->_Finalize();
	  // MOD-BY-LEETEN 01/31/2013-END
	}
};

