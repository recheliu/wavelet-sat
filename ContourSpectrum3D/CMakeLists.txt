# Reference:

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(ContourSpectrum3D)

# The source codes
ADD_EXECUTABLE(${PROJECT_NAME}
	ContourSpectrum3D.cpp
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
	rt	
    ITKNrrdIO	
    z
    )

  if( WITH_BOOST )
    target_link_libraries(${PROJECT_NAME}
      boost_system
      boost_filesystem
      )
  endif()
endif()

if( WITH_NETCDF )
	if( DEFINED WIN32 )
	else()
		target_link_libraries(${PROJECT_NAME}
			hdf5
			netcdf
			)
	endif()
endif()

############################################################
# Copyright (c) 2013 Teng-Yok Lee
#
# See the file LICENSE.txt for copying permission.
############################################################
