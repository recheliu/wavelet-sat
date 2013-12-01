# Reference:
# https://gforge.sci.utah.edu/gf/project/findcuda/scmsvn/?action=browse&path=%2F*checkout*%2Ftrunk%2FFindCuda.html

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(WaveletSAT)

# ADD-BY-LEETEN 09/02/2011-BEGIN
find_package(ThirdPartyLib)

set(MYLIB_DIR 		"~/mylib" 	CACHE PATH "The path to MyLib.")	

# ADD-BY-LEETEN 03/29/2013-BEGIN
set(WITH_DOUBLE_COEF	ON		CACHE BOOL "Is the coefficients stored in double precision?")	
if( WITH_DOUBLE_COEF )
	add_definitions(-DWITH_DOUBLE_COEF=1)
endif()
# ADD-BY-LEETEN 03/29/2013-END
include_directories(
	${LIB_DIR}/include
	${MYLIB_DIR}/include
)

link_directories(
	${LIB_DIR}/lib
	${MYLIB_DIR}/lib
)

# ADD-BY-LEETEN 01/31/2013-BEGIN
if( NOT DEFINED WIN32 )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -O3")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++0x -O3")
endif()
# ADD-BY-LEETEN 01/31/2013-END

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

endif()
# ADD-BY-LEETEN 12/28/2012-END

# ADD-BY-LEETEN 10/28/2012-BEGIN
set(WITH_NETCDF	ON	CACHE 	BOOL "Is NetCDF used?"	FORCE)	
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

# ADD-BY-LEETEN 03/28/2013-BEGIN
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
# ADD-BY-LEETEN 03/28/2013-END

# ADD-BY-LEETEN 01/10/2012-BEGIN
# ADD-BY-LEETEN 01/18/2012-BEGIN
set(WITH_SAT_FILE	OFF	CACHE 	BOOL "Is SAT file format used?")	
if(WITH_SAT_FILE)
	add_definitions(
		-DWITH_SAT_FILE
	)
else()
# ADD-BY-LEETEN 01/18/2012-END

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

	# ADD-BY-LEETEN 2013/12/01-BEGIN

	set(WITHOUT_BIN_AGGREGATION	OFF	CACHE 	BOOL "Is bin aggregation not supported?")	
	if( WITHOUT_BIN_AGGREGATION	) 
		add_definitions(-DWITHOUT_BIN_AGGREGATION)
	endif()	

	# Only build the python binding if the SAT is not written to file.
	add_subdirectory(python)		
	# ADD-BY-LEETEN 2013/12/01-END
# ADD-BY-LEETEN 01/10/2012-END
endif()	# ADD-BY-LEETEN 01/18/2012

add_subdirectory(WaveletSAT)
add_subdirectory(SimpleNDEncoder)
add_subdirectory(ContourSpectrum3D) # ADD-BY-LEETEN 04/20/2013
add_subdirectory(SimpleNDQuery)

add_subdirectory(TestDatatypePerf)	# ADD-BY-LEETEN 10/21/2012

add_subdirectory(Ocean)		# ADD-BY-LEETEN 01/18/2012
add_subdirectory(MJO)		# ADD-BY-LEETEN 03/28/2013

# ADD-BY-LEETEN 01/10/2012-BEGIN
if( WITH_CUDA )
	add_dependencies(SimpleNDEncoder 	CudaDWT)
	
	add_dependencies(Ocean 	CudaDWT)	# ADD-BY-LEETEN 01/18/2012
	add_dependencies(MJO	CudaDWT)	# ADD-BY-LEETEN 03/28/2013
endif()
# ADD-BY-LEETEN 01/10/2012-END

if(NOT WITH_SAT_FILE)	# ADD-BY-LEETEN 03/28/2013-BEGIN
	# ADD-BY-LEETEN 2013/08/03-BEGIN
	set(WITH_VIEWER	OFF	CACHE 	BOOL "Should the viewer be built?")	
	if( WITH_VIEWER )
	# ADD-BY-LEETEN 2013/08/03-END
	add_subdirectory(SATSepDWTViewer2)	
	endif()	# ADD-BY-LEETEN 2013/08/03
endif()					# ADD-BY-LEETEN 03/28/2013-END

# ADD-BY-LEETEN 2013/07/06-BEGIN
configure_file (
  "${PROJECT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  @ONLY
  )
  
install(
	FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
	DESTINATION ".")
# ADD-BY-LEETEN 2013/07/06-END

