# Reference:

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(TestDatatypePerf)

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
	TestDatatypePerf.cpp
)

# ADD-BY-LEETEN 2013/0811-BEGIN
target_link_libraries(${PROJECT_NAME}
	${ThirdPartyLib_LIBRARIES}
	)
# ADD-BY-LEETEN 2013/0811-END

if( DEFINED WIN32 ) 
else()
  target_link_libraries(${PROJECT_NAME}
    rt
    )
endif()

############################################################
# Copyright (c) 2013 Teng-Yok Lee
#
# See the file LICENSE.txt for copying permission.
############################################################
