# Reference:

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(SimpleNDQuery)

find_package(ThirdPartyLib)
set(MYLIB_DIR 		"~/mylib" 	CACHE PATH "The path to MyLib.")	

include_directories(
	${MYLIB_DIR}/include
)

link_directories(
	${MYLIB_DIR}/lib
)

# The source codes
ADD_EXECUTABLE(${PROJECT_NAME}
	SimpleNDQuery.cpp
)

# ADD-BY-LEETEN 2013/0811-BEGIN
target_link_libraries(${PROJECT_NAME}
	${ThirdPartyLib_LIBRARIES}
	)
# ADD-BY-LEETEN 2013/0811-END

if( DEFINED WIN32 ) 
  target_link_libraries(${PROJECT_NAME}
    debug	libopt_d	optimized	libopt_r
    )
else()
  target_link_libraries(${PROJECT_NAME}
    libopt
	rt	# ADD-BY-LEETEN 10/29/2012 
    )
  # ADD-BY-LEETEN 12/31/2012-BEGIN
  if( WITH_BOOST )
    target_link_libraries(${PROJECT_NAME}
      boost_system
      boost_filesystem
      )
  endif()
  # ADD-BY-LEETEN 12/31/2012-END
endif()

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
