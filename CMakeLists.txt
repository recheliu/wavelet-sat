# Reference:
# https://gforge.sci.utah.edu/gf/project/findcuda/scmsvn/?action=browse&path=%2F*checkout*%2Ftrunk%2FFindCuda.html

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(WaveletSAT)

# ADD-BY-LEETEN 09/02/2011-BEGIN
set(LIB_DIR 		"~/lib"		CACHE PATH "The path to 3rd party libraries.")	
set(MYLIB_DIR 		"~/mylib" 	CACHE PATH "The path to MyLib.")	

include_directories(
	${LIB_DIR}/include
	${MYLIB_DIR}/include
)

link_directories(
	${LIB_DIR}/lib
	${MYLIB_DIR}/lib
)

# ADD-BY-LEETEN 12/28/2012-BEGIN
set(WITH_BOOST	OFF	CACHE 	BOOL "Is Boost used?")	
set(BOOST_INC_DIR 	"/usr/local/include" 	CACHE PATH "The path to Boost headers.")	
set(BOOST_LIB_DIR 	"/usr/local/lib" 		CACHE PATH "The path to Boost libraries.")	
if( WITH_BOOST )
	add_definitions(
		-DWITH_BOOST
	)
	include_directories(
		${BOOST_INC_DIR}
	)

	link_directories(
		${BOOST_LIB_DIR}
	)

	# ADD-BY-LEETEN 12/30/2012-BEGIN
	set(WITH_SMART_PTR		OFF	CACHE 	BOOL "Is smart ptr used?")	
	if( WITH_SMART_PTR )
		add_definitions(
			-DWITH_SMART_PTR=1
		)
	endif()
	# ADD-BY-LEETEN 12/30/2012-END
endif()
# ADD-BY-LEETEN 12/28/2012-END

# ADD-BY-LEETEN 10/28/2012-BEGIN
# MOD-BY-LEETEN 01/02/2013-FROM:	set(WITH_NETCDF	ON	CACHE 	BOOL "Is NetCDF used?")	
set(WITH_NETCDF	ON	CACHE 	BOOL "Is NetCDF used?"	FORCE)	
# MOD-BY-LEETEN 01/02/2013-END
set(WITH_NETCDF4	ON	CACHE 	BOOL "Is NetCDF4 used?")	# ADD-BY-LEETEN 11/09/2012
set(NETCDF_DIR "/usr/local" CACHE PATH "The intalled path of NetCDF")	
set(HDF5_DIR "/usr/local" CACHE PATH "The intalled path of HDF5") # ADD-BY-LEETEN 10/29/2012
if( WITH_NETCDF )
	add_definitions(-DWITH_NETCDF)
	# ADD-BY-LEETEN 10/29/2012-BEGIN
	include_directories(
	  "${NETCDF_DIR}/include"
	  "${HDF5_DIR}/include"
	  )
	link_directories(
	  "${NETCDF_DIR}/lib"
	  "${HDF5_DIR}/lib"
	  )
	# ADD-BY-LEETEN 10/29/2012-END
	# ADD-BY-LEETEN 11/09/2012-BEGIN
	if( WITH_NETCDF4 )
	  add_definitions(-DWITH_NETCDF4)
	endif()
	# ADD-BY-LEETEN 11/09/2012-END
endif()
# ADD-BY-LEETEN 10/28/2012-END

include_directories (
	"WaveletSAT"
	"${MYLIB_DIR}/include" )

link_directories (
	"${MYLIB_DIR}/lib" )

# ADD-BY-LEETEN 01/10/2012-BEGIN
set(WITH_CUDA	OFF	CACHE 	BOOL "Is CUDA used?")	
if( WITH_CUDA )
	add_definitions(-DWITH_CUDA)
	
	set(CUDPP_DIR	"/usr/local"	CACHE 	PATH	"The installed path of CUDPP.")	

	FIND_PACKAGE(CUDA)
	include_directories(
		${CUDA_TOOLKIT_ROOT_DIR}/include
		${CUDPP_DIR}/include
		)
	if(CMAKE_SIZEOF_VOID_P EQUAL 8)
		add_definitions(-DWIN64)
		link_directories(
			${CUDA_TOOLKIT_ROOT_DIR}/lib/x64)
	else()
		link_directories(
			${CUDA_TOOLKIT_ROOT_DIR}/lib/Win32)
    endif()	
	link_directories(
		${CUDPP_DIR}/lib
		)

	add_subdirectory(CudaDWT)
	include_directories(
		CudaDWT
		)
		
	link_directories(
		${PROJECT_BINARY_DIR}/CudaDWT
		)
endif()
# ADD-BY-LEETEN 01/10/2012-END

add_subdirectory(WaveletSAT)
add_subdirectory(SimpleND)
# DEL-BY-LEETEN 01/10/2012:	add_dependencies(SimpleND WaveletSAT)

# ADD-BY-LEETEN 09/30/2012-BEGIN
add_subdirectory(SimpleNDWithNrrd)
# DEL-BY-LEETEN 01/10/2012:	add_dependencies(SimpleNDWithNrrd WaveletSAT)
# ADD-BY-LEETEN 09/30/2012-END

# ADD-BY-LEETEN 12/16/2012-BEGIN
add_subdirectory(SimpleNDEncoder)
# DEL-BY-LEETEN 01/10/2012:	add_dependencies(SimpleNDEncoder WaveletSAT)
add_subdirectory(SimpleNDQuery)
# DEL-BY-LEETEN 01/10/2012:	add_dependencies(SimpleNDQuery WaveletSAT)
# ADD-BY-LEETEN 12/16/2012-END

add_subdirectory(TestDatatypePerf)	# ADD-BY-LEETEN 10/21/2012

# ADD-BY-LEETEN 01/10/2012-BEGIN
if( WITH_CUDA )
	add_dependencies(SimpleND 			CudaDWT)
	add_dependencies(SimpleNDWithNrrd 	CudaDWT)
	add_dependencies(SimpleNDEncoder 	CudaDWT)
endif()
# ADD-BY-LEETEN 01/10/2012-END