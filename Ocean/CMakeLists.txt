# Reference:
# This project uses WaveletSAT to analyze the vector field in the ocean dataset
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(Ocean)

set(LIB_DIR 	"~/lib" CACHE PATH "The installed path to 3rd party ibraries.")
include_directories(${LIB_DIR}/include)
link_directories(${LIB_DIR}/lib)

set(HDF5_DIR	"~/lib" CACHE PATH "The installed path to HDF5.")
include_directories(${HDF5_DIR}/include)
link_directories(${HDF5_DIR}/lib)

set(NETCDF_DIR 	"~/lib" CACHE PATH "The installed path to NetCDF 4.")
include_directories(${NETCDF_DIR}/include)
link_directories(${NETCDF_DIR}/lib)

set(MYLIB_DIR 		"~/mylib" CACHE PATH "The installed path to my own libraries.")
include_directories(${MYLIB_DIR}/include)
link_directories(${MYLIB_DIR}/lib)

# set(WAVELETSAT_DIR 	"~/mylib" CACHE PATH "The installed path of WaveletSAT.")
# include_directories(${WAVELETSAT_DIR}/include)
# link_directories(${WAVELETSAT_DIR}/lib)

set(WITH_SAT_FILE	OFF	CACHE BOOL	"Whether SAT is used.")
if( WITH_SAT_FILE )
	add_definitions(-DWITH_SAT_FILE)
else()
	set(WITH_CUDA 	OFF	CACHE BOOL	"Whether CUDA is used or not.")
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
	endif()
endif()

# add_definitions(
	# -DWITH_NETCDF
	# -DWITH_NETCDF4
# )

ADD_EXECUTABLE(${PROJECT_NAME}
	Vector2D.h
	WaveletSAT_Ocean_main.cpp
)

target_link_libraries(
	${PROJECT_NAME}
	debug	libopt_d	optimized	libopt_r
	)

if( NOT WITH_SAT_FILE AND WITH_CUDA )
	target_link_libraries(
		${PROJECT_NAME}
		debug	CudaDWT_d	optimized	CudaDWT_r
		)
endif()