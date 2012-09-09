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
    libopt
    )
endif()
# ADD-BY-LEETEN 09/09/2012-END
# ADD-BY-LEETEN 09/09/2012-END
