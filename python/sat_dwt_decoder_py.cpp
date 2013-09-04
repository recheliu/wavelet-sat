#include <windows.h>
#include <iostream>
#include <list>
#include <vector>
using namespace std;
#include <boost/python.hpp>
#include "boost/python/class.hpp" 
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/args.hpp>	// ADD-BY-LEETEN 2013/09/03
#include <boost/python/extract.hpp>
#include <boost/python/to_python_value.hpp>
using namespace boost::python;

#include "libclock.h"

#include "SimpleNDFile.h"
using namespace WaveletSAT;

// Ref: http://wiki.python.org/moin/boost.python/PointersAndSmartPointers

// MOD-BY-LEETEN 2013/09/03-FROM:	struct pySATSepDWTDecoder: public CSimpleNDFile<double, double, WaveletSAT::typeBin, double> {
struct pySATSepDWTDecoder: public CSimpleNDFile<double, WaveletSAT::typeSum, WaveletSAT::typeBin, WaveletSAT::typeWavelet> {
// MOD-BY-LEETEN 2013/09/03-END

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

	template<class T>
	void 
	_ConvertVectorToList
	(
		const vector<T>& vVector,
		boost::python::list& python_list
	)
	{
		while( boost::python::len(python_list) > 0 )
			python_list.pop();

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

	#if	0	// MOD-BY-LEETEN 2013/09/03-FROM:
	boost::python::list
	get_size
	(
	) 
	{
		vector<size_t> vuDataSize;
		this->_GetDataSize(vuDataSize);
		boost::python::list* data_size = new boost::python::list();
		for(size_t d = 0; d < this->UGetNrOfDims(); d++) 
			data_size->append(vuDataSize[d]);
		return *data_size;
	}
	#else	// MOD-BY-LEETEN 2013/09/03-TO:
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
	#endif	// MOD-BY-LEETEN 2013/09/03-END
	
	void	
	load
	(
		std::string filepath
	)		
	{	
		_LoadFile(filepath.c_str());	
	}

	#if	0	// MOD-BY-LEETEN 2013/09/03-FROM:
	boost::python::list
	get_region_histogram
	(
		boost::python::list region_left, 
		boost::python::list region_right
	)		
	#else	// MOD-BY-LEETEN 2013/09/03-TO:
	void
	get_region_histogram
	(
		const boost::python::list& region_left, 
		const boost::python::list& region_right,
		boost::python::list& region_hist
	)
	#endif	// MOD-BY-LEETEN 2013/09/03-END
	{	
		#if	0	// MOD-BY-LEETEN 2013/09/03-FROM:
		vector<size_t> vuLeft;
		vuLeft.resize((size_t)boost::python::len(region_left));
		for(size_t d = 0; d < vuLeft.size(); d++)
			vuLeft[d] = boost::python::extract<size_t>(region_left[d]);
			
		vector<size_t> vuRight;
		vuRight.resize((size_t)boost::python::len(region_right));
		for(size_t d = 0; d < vuRight.size(); d++)
			vuRight[d] = boost::python::extract<size_t>(region_right[d]);
		#else	// MOD-BY-LEETEN 2013/09/03-TO:
		vector<size_t> vuLeft;
		_ConvertListToVector(region_left, vuLeft);

		vector<size_t> vuRight;
		_ConvertListToVector(region_right, vuRight);
		#endif	// MOD-BY-LEETEN 2013/09/03-END
		vector<WaveletSAT::typeSum> vdHist;
		_GetRegionSums(vuLeft, vuRight, vdHist);
		
		#if	0	// MOD-BY-LEETEN 2013/09/03-FROM:
		boost::python::list* region_hist = new boost::python::list();
		for(size_t b = 0; b < vdHist.size(); b++)
			region_hist->append(vdHist[b]);

		return *region_hist;
		#else	// MOD-BY-LEETEN 2013/09/03-TO:
		_ConvertVectorToList<WaveletSAT::typeSum>(vdHist, region_hist);
		#endif	// MOD-BY-LEETEN 2013/09/03-END
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

		.def("get_region_histogram", &pySATSepDWTDecoder::get_region_histogram)
		;	
}



