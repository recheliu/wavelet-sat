# Reference:
# https://gforge.sci.utah.edu/gf/project/findcuda/scmsvn/?action=browse&path=%2F*checkout*%2Ftrunk%2FFindCuda.html

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(WaveletSAT)

find_package(ThirdPartyLib)

set(MYLIB_DIR 		"~/mylib" 	CACHE PATH "The path to MyLib.")	

set(WITH_DOUBLE_COEF	ON		CACHE BOOL "Is the coefficients stored in double precision?")	
if( WITH_DOUBLE_COEF )
	add_definitions(-DWITH_DOUBLE_COEF=1)
endif()
include_directories(
	${LIB_DIR}/include
	${MYLIB_DIR}/include
)

link_directories(
	${LIB_DIR}/lib
	${MYLIB_DIR}/lib
)

if( NOT DEFINED WIN32 )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -O3")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++0x -O3")
endif()

# # DEL-BY-LEETEN 2015/04/30-BEGIN
# set(WITH_BOOST	OFF	CACHE 	BOOL "Is Boost used?")	
# set(BOOST_INC_DIR 	"/usr/local/include" 	CACHE PATH "The path to Boost headers.")	
# set(BOOST_LIB_DIR 	"/usr/local/lib" 		CACHE PATH "The path to Boost libraries.")	
# if( WITH_BOOST )
	# add_definitions(
		# -DWITH_BOOST
	# )
	# include_directories(
		# ${BOOST_INC_DIR}
	# )

	# link_directories(
		# ${BOOST_LIB_DIR}
	# )

# endif()
# # DEL-BY-LEETEN 2015/04/30-END

set(WITH_NETCDF	ON	CACHE 	BOOL "Is NetCDF used?"	FORCE)	
set(WITH_NETCDF4	ON	CACHE 	BOOL "Is NetCDF4 used?")	
set(NETCDF_DIR "/usr/local" CACHE PATH "The intalled path of NetCDF")	
set(HDF5_DIR "/usr/local" CACHE PATH "The intalled path of HDF5") 
if( WITH_NETCDF )
	add_definitions(-DWITH_NETCDF)
	include_directories(
	  "${NETCDF_DIR}/include"
	  "${HDF5_DIR}/include"
	  )
	link_directories(
	  "${NETCDF_DIR}/lib"
	  "${HDF5_DIR}/lib"
	  )
	if( WITH_NETCDF4 )
	  add_definitions(-DWITH_NETCDF4)
	endif()
endif()

include_directories (
	"WaveletSAT"
	"${MYLIB_DIR}/include" )

link_directories (
	"${MYLIB_DIR}/lib" )

FIND_PACKAGE(CUDA)
include_directories(
	${CUDA_TOOLKIT_ROOT_DIR}/include
	)
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
	add_definitions(-DWIN64)
	link_directories(
		${CUDA_TOOLKIT_ROOT_DIR}/lib/x64)
else()
	link_directories(
		${CUDA_TOOLKIT_ROOT_DIR}/lib/Win32)
endif()	

set(WITH_SAT_FILE	OFF	CACHE 	BOOL "Is SAT file format used?")	
if(WITH_SAT_FILE)
	add_definitions(
		-DWITH_SAT_FILE
	)
else()

set(WITH_CUDA	OFF	CACHE 	BOOL "Is CUDA used?")	
if( WITH_CUDA )
	add_definitions(-DWITH_CUDA)
	
	add_subdirectory(CudaDWT)
	include_directories(
		CudaDWT
		)
		
	link_directories(
		${PROJECT_BINARY_DIR}/CudaDWT
		)
endif()

	set(WITHOUT_BIN_AGGREGATION	OFF	CACHE 	BOOL "Is bin aggregation not supported?")	
	if( WITHOUT_BIN_AGGREGATION	) 
		add_definitions(-DWITHOUT_BIN_AGGREGATION)
	endif()	

	# Only build the python binding if the SAT is not written to file.
	# MOD-BY-LEETEN 2015/04/29:	add_subdirectory(python)		
	# TODO: Use numpy to re-implement the python binding.
	# MOD-BY-LEETEN 2015/04/29-END
	
	# ADD-BY-LEETEN 2015/04/30-BEGIN
	set(WITH_PYTHON	OFF	CACHE 	BOOL "Building the python binding?")	
	IF(WITH_PYTHON)
		add_subdirectory(python)		
	ENDIF()
	# ADD-BY-LEETEN 2015/04/30-END
endif()	

add_subdirectory(WaveletSAT)
add_subdirectory(SimpleNDEncoder)
add_subdirectory(ContourSpectrum3D) 
add_subdirectory(SimpleNDQuery)

add_subdirectory(TestDatatypePerf)	

add_subdirectory(Ocean)		
add_subdirectory(MJO)		

if( WITH_CUDA )
	add_dependencies(SimpleNDEncoder 	CudaDWT)
	
	add_dependencies(Ocean 	CudaDWT)	
	add_dependencies(MJO	CudaDWT)	
endif()

if(NOT WITH_SAT_FILE)	
	set(WITH_VIEWER	OFF	CACHE 	BOOL "Should the viewer be built?")	
	if( WITH_VIEWER )
		add_subdirectory(SATSepDWTViewer2)	
	endif()	
endif()

configure_file (
  "${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  @ONLY
  )
  
install(
	FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
	DESTINATION ".")


############################################################
# Copyright (c) 2013 Teng-Yok Lee
#
# See the file LICENSE.txt for copying permission.
############################################################
