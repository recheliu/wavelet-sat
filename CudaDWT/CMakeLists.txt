# Reference:
# https://gforge.sci.utah.edu/gf/project/findcuda/scmsvn/?action=browse&path=%2F*checkout*%2Ftrunk%2FFindCuda.html

#The name of the project
PROJECT(CudaDWT)

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#Searching CUDA
FIND_PACKAGE(CUDA)

CUDA_INCLUDE_DIRECTORIES(
	${CUDA_CUT_INCLUDE_DIRS} 
	${CUDA_INCLUDE_DIRS}
	)

# The source codes
CUDA_ADD_LIBRARY(${PROJECT_NAME}
	CudaDWT.h
	CudaDWT.cu
	ProjToWavelet_kernel.h
	MarkSegments_kernel.h
)

set_target_properties(${PROJECT_NAME} PROPERTIES 
	DEBUG_OUTPUT_NAME "${PROJECT_NAME}_d"
	RELEASE_OUTPUT_NAME "${PROJECT_NAME}_r"
	)

# ADD-BY-LEETEN 01/18/2012-BEGIN
install(
	TARGETS ${PROJECT_NAME}
	ARCHIVE	DESTINATION lib
	)

install(
	FILES CudaDWT.h
	DESTINATION include
	)
# ADD-BY-LEETEN 01/18/2012-END