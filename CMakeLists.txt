# Reference:
# https://gforge.sci.utah.edu/gf/project/findcuda/scmsvn/?action=browse&path=%2F*checkout*%2Ftrunk%2FFindCuda.html

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(WaveletSAT)

# ADD-BY-LEETEN 09/02/2011-BEGIN
# MOD-BY-LEETEN 12/22/2012-FROM:
	#	set(MYLIB_DIR "~/mylib" CACHE PATH "The path to my own libraries")
# TO:
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
# MOD-BY-LEETEN 12/22/2012-END

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
set(WITH_NETCDF	ON	CACHE 	BOOL "Is NetCDF used?")	
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
	# MOD-BY-LEETEN 10/27/2012-FROM:	"SAT"
	"WaveletSAT"
	# MOD-BY-LEETEN 10/27/2012-END
	"${MYLIB_DIR}/include" )

link_directories (
	"${MYLIB_DIR}/lib" )

# DEL-BY-LEETEN 12/16/2012:	add_subdirectory(SAT)
add_subdirectory(WaveletSAT)
add_subdirectory(SimpleND)
# MOD-BY-LEETEN 10/27/2012-FROM:	add_dependencies(SimpleND SAT)
add_dependencies(SimpleND WaveletSAT)
# MOD-BY-LEETEN 10/27/2012-END

# ADD-BY-LEETEN 09/30/2012-BEGIN
add_subdirectory(SimpleNDWithNrrd)
# MOD-BY-LEETEN 10/27/2012-FROM:	add_dependencies(SimpleNDWithNrrd SAT)
add_dependencies(SimpleNDWithNrrd WaveletSAT)
# MOD-BY-LEETEN 10/27/2012-END
# ADD-BY-LEETEN 09/30/2012-END

# ADD-BY-LEETEN 12/16/2012-BEGIN
add_subdirectory(SimpleNDEncoder)
add_dependencies(SimpleNDEncoder WaveletSAT)
add_subdirectory(SimpleNDQuery)
add_dependencies(SimpleNDQuery WaveletSAT)
# ADD-BY-LEETEN 12/16/2012-END

add_subdirectory(TestDatatypePerf)	# ADD-BY-LEETEN 10/21/2012

