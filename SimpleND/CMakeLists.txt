# Reference:

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(SimpleND)

set(MYLIB_DIR 		"~/mylib" 	CACHE PATH "The path to MyLib.")	

include_directories(
	${MYLIB_DIR}/include
)

link_directories(
	${MYLIB_DIR}/lib
)

# The source codes
ADD_EXECUTABLE(${PROJECT_NAME}
	SimpleND.cpp
)

# ADD-BY-LEETEN 09/09/2012-BEGIN
if( DEFINED WIN32 ) # ADD-BY-LEETEN 09/09/2012
  target_link_libraries(${PROJECT_NAME}
    debug	libopt_d	optimized	libopt_r
    )
# ADD-BY-LEETEN 09/09/2012-BEGIN
else()
  target_link_libraries(${PROJECT_NAME}
	rt	# ADD-BY-LEETEN 10/29/2012 
    libopt
    )
  # ADD-BY-LEETEN 01/02/2013-BEGIN
  if( WITH_BOOST )
    target_link_libraries(${PROJECT_NAME}
      boost_system
      boost_filesystem
      )
  endif()
  # ADD-BY-LEETEN 01/02/2013-END
endif()
# ADD-BY-LEETEN 09/09/2012-END
# ADD-BY-LEETEN 09/09/2012-END

# ADD-BY-LEETEN 10/28/2012-BEGIN
if( WITH_NETCDF )
	if( DEFINED WIN32 )
	else()
		target_link_libraries(${PROJECT_NAME}
			hdf5
			netcdf
			)
	endif()
endif()
# ADD-BY-LEETEN 10/28/2012-END
