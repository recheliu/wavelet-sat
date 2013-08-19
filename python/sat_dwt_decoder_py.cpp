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
#include <boost/python/extract.hpp>
#include <boost/python/to_python_value.hpp>
using namespace boost::python;

#include "libclock.h"

#include "SimpleNDFile.h"
using namespace WaveletSAT;

// Ref: http://wiki.python.org/moin/boost.python/PointersAndSmartPointers

struct pySATSepDWTDecoder: public CSimpleNDFile<double, double, WaveletSAT::typeBin, double> {
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
	
	void	
	load
	(
		std::string filepath
	)		
	{	
		_LoadFile(filepath.c_str());	
	}

	boost::python::list
	get_region_histogram
	(
		boost::python::list region_left, 
		boost::python::list region_right
	)		
	{	
		vector<size_t> vuLeft;
		vuLeft.resize((size_t)boost::python::len(region_left));
		for(size_t d = 0; d < vuLeft.size(); d++)
			vuLeft[d] = boost::python::extract<size_t>(region_left[d]);
			
		vector<size_t> vuRight;
		vuRight.resize((size_t)boost::python::len(region_right));
		for(size_t d = 0; d < vuRight.size(); d++)
			vuRight[d] = boost::python::extract<size_t>(region_right[d]);

		vector<WaveletSAT::typeSum> vdHist;
		_GetRegionSums(vuLeft, vuRight, vdHist);
		
		boost::python::list* region_hist = new boost::python::list();
		for(size_t b = 0; b < vdHist.size(); b++)
			region_hist->append(vdHist[b]);

		return *region_hist;
	}
};

BOOST_PYTHON_MODULE(sat_dwt_decoder)
{
	//        .def("create",&pySATSepDWTDecoder::create )
	// A method to get the histograms 
    class_<pySATSepDWTDecoder>("decoder")
		.def("load", &pySATSepDWTDecoder::load)
		.def("get_size", &pySATSepDWTDecoder::get_size)
		.def("get_region_histogram", &pySATSepDWTDecoder::get_region_histogram)
		;	
}



