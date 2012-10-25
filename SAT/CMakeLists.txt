# Reference:

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
# MOD-BY-LEETEN 10/17/2012-FROM:
	# PROJECT(WaveletSATLib)
# TO:
PROJECT(SAT)
# MOD-BY-LEETEN 10/17/2012-END

set(MYLIB_DIR 		"~/mylib" 	CACHE PATH "The path to MyLib.")	

include_directories(
	${MYLIB_DIR}/include
)

link_directories(
	${MYLIB_DIR}/lib
)

# The source codes
set(HEADER_FILES
		WaveletSAT.h
		Base.h	# ADD-BY-LEETEN 09/29/2012
		# ADD-BY-LEETEN 09/30/2012-BEGIN
		IsotropicWaveletSAT.h	
		SimpleND.h	
		# ADD-BY-LEETEN 09/30/2012-END
		IntegralHistogram.h			# ADD-BY-LEETEN 10/17/2012
		IntegralHistogramSDWT		# ADD-BY-LEETEN 10/25/2012
	)

add_custom_target (
	${PROJECT_NAME}
	SOURCES
		${HEADER_FILES}
)

install(
	FILES ${HEADER_FILES}
	DESTINATION "include")
	