# Reference:
# https://gforge.sci.utah.edu/gf/project/findcuda/scmsvn/?action=browse&path=%2F*checkout*%2Ftrunk%2FFindCuda.html

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(WaveletSAT)

# ADD-BY-LEETEN 09/02/2011-BEGIN
set(MYLIB_DIR "~/mylib" CACHE PATH "The path to my own libraries")

include_directories (
	"SAT"
	"${MYLIB_DIR}/include" )

link_directories (
	"${MYLIB_DIR}/lib" )

add_subdirectory(SAT)
add_subdirectory(SimpleND)

add_dependencies(SimpleND SAT)

# ADD-BY-LEETEN 09/30/2012-BEGIN
add_subdirectory(SimpleNDWithNrrd)
add_dependencies(SimpleNDWithNrrd SAT)
# ADD-BY-LEETEN 09/30/2012-END


