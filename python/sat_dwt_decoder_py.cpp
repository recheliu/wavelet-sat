#include <windows.h>
#include <iostream>
#include <list>
#include <vector>
using namespace std;
#include <Python.h>	// ADD-BY-LEETEN 2013/09/16
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

		.def("get_region_histogram", &pySATSepDWTDecoder::get_region_histogram)
		;	
}



