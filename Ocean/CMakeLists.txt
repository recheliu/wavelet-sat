# Reference:
# This project uses WaveletSAT to analyze the vector field in the ocean dataset
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(Ocean)

find_package(ThirdPartyLib)

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

# ADD-BY-LEETEN 2013/0811-BEGIN
target_link_libraries(${PROJECT_NAME}
	${ThirdPartyLib_LIBRARIES}
	)
# ADD-BY-LEETEN 2013/0811-END
if( DEFINED WIN32) # ADD-BY-LEETEN 01/31/2013
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
# ADD-BY-LEETEN 01/31/2013-BEGIN
else()
target_link_libraries(
	${PROJECT_NAME}
	rt
	boost_system
	netcdf
	libopt
	)

      if( NOT WITH_SAT_FILE AND WITH_CUDA )
	target_link_libraries(
		${PROJECT_NAME}
		CudaDWT
		)
	    endif()
endif()
# ADD-BY-LEETEN 01/31/2013-END
