# Reference:

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(TestDatatypePerf)

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
	TestDatatypePerf.cpp
)

if( DEFINED WIN32 ) 
else()
  target_link_libraries(${PROJECT_NAME}
    rt
    )
endif()
