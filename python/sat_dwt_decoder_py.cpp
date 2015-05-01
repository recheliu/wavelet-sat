#include <windows.h>
#include <iostream>
#include <list>
#include <vector>
using namespace std;
// ADD-BY-LEETEN 2015/04/30-BEGIN
// Do not enable debug mode for Python 
// since python27_d.lib might be unavailable.
#if	defined(_DEBUG)
	#define WITH_DEBUG	
	#undef _DEBUG
#endif
// ADD-BY-LEETEN 2015/04/30-END
#include <Python.h>	// ADD-BY-LEETEN 2013/09/16
// ADD-BY-LEETEN 2015/04/30-BEGIN
#if	defined(WITH_DEBUG)
	#define _DEBUG
#endif
// ADD-BY-LEETEN 2015/04/30-END
#include <boost/python.hpp>
#include "boost/python/class.hpp" 
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/args.hpp>	// ADD-BY-LEETEN 2013/09/03
#include <boost/python/call.hpp>	// ADD-BY-LEETEN 2013/09/16
#include <boost/python/extract.hpp>
#include <boost/python/to_python_value.hpp>
using namespace boost::python;

#include "libclock.h"

#include "SimpleNDFile.h"
using namespace WaveletSAT;

// Ref: http://wiki.python.org/moin/boost.python/PointersAndSmartPointers

// ADD-BY-LEETEN 2014/01/08-BEGIN
struct pyVector: public vector<double> {
	double 
	get(
		const size_t& index
	)
	{
		const vector<double>& vec = *this;
		return vec[index];
	}

	void 
	set(
		const size_t& index,
		double value
	)
	{
		vector<double>& vec = *this;
		vec[index] = value;
	}
};
// ADD-BY-LEETEN 2014/01/08-END

struct pySATSepDWTDecoder: public CSimpleNDFile<double, WaveletSAT::typeSum, WaveletSAT::typeBin, WaveletSAT::typeWavelet> {

	// ADD-BY-LEETEN 2013/09/03-BEGIN
	template<class T>
	void 
	_ConvertListToVector
	(
		const boost::python::list& python_list, 
		vector<T>& vVector
	)
	{
		size_t uLen = boost::python::len(python_list);
		vVector.resize(uLen);
		for(size_t d = 0; d < uLen; d++)
			vVector[d] = boost::python::extract<T>(python_list[d]);
	}

	// ADD-BY-LEETEN 2013/09/16-BEGIN
	void 
	_ClearList
	(
		boost::python::list& python_list
	)
	{
		while( boost::python::len(python_list) > 0 )
			python_list.pop();
	}
	// ADD-BY-LEETEN 2013/09/16-END

	template<class T>
	void 
	_ConvertVectorToList
	(
		const vector<T>& vVector,
		boost::python::list& python_list
	)
	{
		// ADD-BY-LEETEN 2013/09/16-BEGIN
		_ClearList(python_list);
		// ADD-BY-LEETEN 2013/09/16-END

		for(size_t d = 0; d < vVector.size(); d++)
			python_list.append(vVector[d]);
	}

	int 
	get_n_bins
	(
	) 
	{
		return (int)this->UGetNrOfBins();
	}
	// ADD-BY-LEETEN 2013/09/03-END

	void
	get_size
	(
		boost::python::list& data_size
	) 
	{
		vector<size_t> vuDataSize = this->VGetDimLengths();
		_ConvertVectorToList<size_t>(vuDataSize, data_size);
	}

	void
	get_level
	(
		boost::python::list& data_level
	) 
	{
		vector<size_t> vuDataLevel = VGetDimLevels();
		_ConvertVectorToList<size_t>(vuDataLevel, data_level);
	}
	
	void	
	load
	(
		std::string filepath
	)		
	{	
		_LoadFile(filepath.c_str());	
	}

	void
	get_region_histogram
	(
		const boost::python::list& region_left, 
		const boost::python::list& region_right,
		boost::python::list& region_hist
	)
	{	
		vector<size_t> vuLeft;
		_ConvertListToVector(region_left, vuLeft);

		vector<size_t> vuRight;
		_ConvertListToVector(region_right, vuRight);
		vector<WaveletSAT::typeSum> vdHist;
		_GetRegionSums(vuLeft, vuRight, vdHist);
		
		_ConvertVectorToList<WaveletSAT::typeSum>(vdHist, region_hist);
	}

	// ADD-BY-LEETEN 2013/09/03-BEGIN
	void get_coef_sparse
	(
		const boost::python::list& level, 
		const boost::python::list& local_coord,
		boost::python::dict& coef
	) 
	{
		vector<size_t> vuLevel;
		_ConvertListToVector<size_t>(level, vuLevel);

		vector<size_t> vuLocalCoord;
		_ConvertListToVector<size_t>(local_coord, vuLocalCoord);

		vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> > vpairCoefBinValues;
		_GetCoefSparse
		(
			vuLevel, 
			vuLocalCoord, 
			vpairCoefBinValues
		);

		coef.clear();
		for(vector< pair<WaveletSAT::typeBin, WaveletSAT::typeWavelet> >::iterator 
				ivpairCoefBinValues = vpairCoefBinValues.begin();
			ivpairCoefBinValues != vpairCoefBinValues.end();
			ivpairCoefBinValues ++) 
		{
			coef[ivpairCoefBinValues->first] = ivpairCoefBinValues->second;
		}
	}
	// ADD-BY-LEETEN 2013/09/03-END

	// ADD-BY-LEETEN 2013/09/16-BEGIN
	void apply_filter
	(
		PyObject* callable,
		const boost::python::list& left_offset, 
		const boost::python::list& right_offset,
		boost::python::list& result,
		const boolean verbose
	) 
	{
		LIBCLOCK_INIT(verbose, __FUNCTION__);

		LIBCLOCK_BEGIN(verbose);
		vector<int> viLeftOffset;
		_ConvertListToVector<int>(left_offset, viLeftOffset);
		vector<size_t> vuLeft;
		vuLeft.resize(this->UGetNrOfDims());

		vector<int> viRightOffset;
		_ConvertListToVector<int>(right_offset, viRightOffset);
		vector<size_t> vuRight;
		vuRight.resize(this->UGetNrOfDims());

		vector<size_t> vuCenter;
		vuCenter.resize(this->UGetNrOfDims());

		vector<double> vdResult;
		vdResult.resize(uDataSize);
		LIBCLOCK_END(verbose);

		LIBCLOCK_BEGIN(verbose);
		for(size_t v = 0; v < this->uDataSize; v++) {
			WaveletSAT::_ConvertIndexToSub(v, vuCenter, this->vuDimLengths);
			for(size_t d = 0; d < this->UGetNrOfDims(); d++) {
				int iLeft = (int)vuCenter[d] + viLeftOffset[d];
				vuLeft[d] = (size_t)min(max(iLeft, 0), (int)vuDimLengths[d] - 1);

				int iRight = (int)vuCenter[d] + viRightOffset[d];
				vuRight[d] = (size_t)min(max(iRight, 0), (int)vuDimLengths[d] - 1);
			}
			vector<WaveletSAT::typeSum> vdHist;
			_GetRegionSums(vuLeft, vuRight, vdHist);

			boost::python::list hist;
			_ConvertVectorToList<WaveletSAT::typeSum>(vdHist, hist);

			vdResult[v] = call<double, boost::python::list>(callable, boost::ref(hist));
		}
		LIBCLOCK_END(verbose);

		LIBCLOCK_BEGIN(verbose);
		_ConvertVectorToList<double>(vdResult, result);
		LIBCLOCK_END(verbose);

		LIBCLOCK_PRINT(verbose);
	}
	// ADD-BY-LEETEN 2013/09/16-END

	#if	1	// TEST-ADD
	/*
	void
	clamp_to_datasize
	(
		const boost::python::list& src,
		boost::python::list& dst,
		const boolean verbose
	)
	{
		vector<double> vdSrc;
		_ConvertListToVector<double>(src, vdSrc);
		vector<double> vdDst;
		_ClampToDataSize(vdSrc, vdDst);
		_ConvertVectorToList<double>(vdDst, dst);
	}

	void
	comp_region_sum
	(
		const boost::python::list& left_offset, 
		const boost::python::list& right_offset,
		const boost::python::list& sat,
		const boost::python::list& sat_lengths,	
		boost::python::list& region_sum,
		const boolean verbose
	)
	{
		vector<int> viLeftOffset;		
		_ConvertListToVector<int>(left_offset, viLeftOffset);
		vector<int> viRightOffset;		
		_ConvertListToVector<int>(right_offset, viRightOffset);
		vector<double> vdSAT;			
		_ConvertListToVector<double>(sat, vdSAT);
		vector<size_t> vuSATLengths;		
		_ConvertListToVector<size_t>(sat_lengths, vuSATLengths);

		vector<double> vdRegionSum;
		vdRegionSum.resize(vdSAT.size());
		_ComputeRegionSum
		(
			viLeftOffset,
			viRightOffset,
			vdSAT,
			vuSATLengths,
			vdRegionSum
		);
		_ConvertVectorToList<double>(vdRegionSum, region_sum);
	}

	void 
	get_bin_aggregate
	(
		const int bin_left,
		const int bin_right,
		boost::python::list& result,
		const boolean verbose
	) 
	{
		vector<double> vdResult;
		vdResult.resize(uDataSize);
		_DecodeAggregatedBin(bin_left, bin_right, vdResult);
		_ConvertVectorToList<double>(vdResult, result);
	}
	*/
	void
	clamp_to_datasize
	(
		const pyVector& src,
		pyVector& dst,
		const boolean verbose
	)
	{
		_ClampToDataSize(src, dst);
	}

	void
	comp_region_sum
	(
		const boost::python::list& left_offset, 
		const boost::python::list& right_offset,
		const pyVector& sat,
		const boost::python::list& sat_lengths,	
		pyVector& region_sum,
		const boolean verbose
	)
	{
		vector<int> viLeftOffset;		
		_ConvertListToVector<int>(left_offset, viLeftOffset);
		vector<int> viRightOffset;		
		_ConvertListToVector<int>(right_offset, viRightOffset);
		vector<size_t> vuSATLengths;		
		_ConvertListToVector<size_t>(sat_lengths, vuSATLengths);

		_ComputeRegionSum
		(
			viLeftOffset,
			viRightOffset,
			sat,
			vuSATLengths,
			region_sum
		);
	}

	void 
	clamp_border
	(
		const boost::python::list& left_offset, 
		const boost::python::list& right_offset,
		pyVector& sat,
		const boolean verbose
	) 
	{
		vector<int> viLeftOffset;		
		_ConvertListToVector<int>(left_offset, viLeftOffset);
		vector<int> viRightOffset;		
		_ConvertListToVector<int>(right_offset, viRightOffset);
		_ClampBorder(sat, viLeftOffset, viRightOffset);
	}

	void 
	get_bin_aggregate
	(
		const int bin_left,
		const int bin_right,
		pyVector& bin_aggregate,
		const boolean verbose
	) 
	{
		_DecodeAggregatedBin(bin_left, bin_right, bin_aggregate);
	}
#endif
};

BOOST_PYTHON_MODULE(sat_dwt_decoder)
{
	//        .def("create",&pySATSepDWTDecoder::create )
	// A method to get the histograms 
    class_<pySATSepDWTDecoder>("decoder")
		.def("load", &pySATSepDWTDecoder::load)
		.def("get_size", &pySATSepDWTDecoder::get_size)

		// ADD-BY-LEETEN 2013/09/03-BEGIN
		.def("get_level", &pySATSepDWTDecoder::get_level,
			args("data_level"), 
			"Get the levels of wavelet coefficients in all dimensions.")

		.def("get_coef_sparse", &pySATSepDWTDecoder::get_coef_sparse, 
			args("level", "local_coord", "sparse_coef"), 
			"Get the wavelet coefficients at the specified level and coordinates within that level.")

		.def("get_n_bins", &pySATSepDWTDecoder::get_n_bins,
			"Get #Bins.")
		// ADD-BY-LEETEN 2013/09/03-END

		// ADD-BY-LEETEN 2013/09/16-BEGIN
		.def("apply_filter", &pySATSepDWTDecoder::apply_filter,
			args("callable", "left_offset", "right_offset", "result", "verbose"), 
			"Apply the filter to all region histograms of the given size.")
		// ADD-BY-LEETEN 2013/09/16-END

#if	1	// TEST-ADD-BEGIN
		.def("comp_region_sum", &pySATSepDWTDecoder::comp_region_sum,
			args("left_offset", "right_offset", "sat", "sat_lengths", "region_sum", "verbose"), 
			"Get the region sum for the given SAT.")

		.def("clamp_to_datasize", &pySATSepDWTDecoder::clamp_to_datasize,
			args("src", "dst", "verbose"), 
			"Clamp the given field to the data size.")

		.def("clamp_border", &pySATSepDWTDecoder::clamp_border,
			args("left_offset", "right_offset", "sat", "verbose"), 
			"Clamp the border according to the left_offset and right_offset.")

		.def("get_bin_aggregate", &pySATSepDWTDecoder::get_bin_aggregate,
			args("bin_left", "bin_right", "result", "verbose"), 
			"Get the sum of bin SATs with bins from the left edge of bin_left to the right edge of bin_right.")
#endif

		.def("get_region_histogram", &pySATSepDWTDecoder::get_region_histogram)
		;	

#if	1	// TEST-ADD
    class_<pyVector>("std_vector")
		.def("get", &pyVector::get,
			args("index"), 
			"Get the value at the given index.")

		.def("set", &pyVector::set,
			args("index", "value"), 
			"Set the value at the given index.")
		;
#endif
}



