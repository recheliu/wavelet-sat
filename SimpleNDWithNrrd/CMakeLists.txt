# Reference:

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(SimpleNDWithNrrd)

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

# The source codes
ADD_EXECUTABLE(${PROJECT_NAME}
	SimpleNDWithNrrd.cpp
)

if( DEFINED WIN32 ) 
  target_link_libraries(${PROJECT_NAME}
    debug	libopt_d	optimized	libopt_r
    debug	ITKNrrdIOd	optimized	ITKNrrdIO
    debug	zlibd		optimized	zlib
    )
else()
  target_link_libraries(${PROJECT_NAME}
    libopt
	rt	# ADD-BY-LEETEN 10/29/2012 
    ITKNrrdIO	
    z
    )
endif()

# ADD-BY-LEETEN 10/28/2012-BEGIN
if( WITH_NETCDF )
	if( DEFINED WIN32 )
		# # DEL-BY-LEETEN 12/22/2012-BEGIN
		# target_link_libraries(${PROJECT_NAME}
			# libszip.lib
			# zlib.lib
			# hdf5.lib
			# hdf5_cpp.lib
			# hdf5_hl.lib
			# hdf5_hl_cpp.lib
			# netcdf.lib
			# )
		# # DEL-BY-LEETEN 12/22/2012-END
	else()
		target_link_libraries(${PROJECT_NAME}
			# DEL-BY-LEETEN 10/28/2012: szip
			# DEL-BY-LEETEN 10/28/2012: zlib
			hdf5
			netcdf
			)
	endif()
endif()
# ADD-BY-LEETEN 10/28/2012-END
