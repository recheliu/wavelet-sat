# Reference:

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(WaveletSAT)

set(WITH_GNUPLOT OFF			CACHE BOOL	"The flag whethre GNUPLOT is enabled.")
set(GNUPLOT_BIN_DIR "/usr/bin"	CACHE PATH	"The path to the executable of GNU Plot.")
if(WITH_GNUPLOT)
	set(IS_USING_GNUPLOT	1)
else()
	set(IS_USING_GNUPLOT	0)
endif()
if(DEFINED WIN32)	
	set(GNUPLOT_BIN_PATHFILENAME	${GNUPLOT_BIN_DIR}/gnuplot.exe)
else()
	set(GNUPLOT_BIN_PATHFILENAME	${GNUPLOT_BIN_DIR}/gnuplot)
endif()
configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/gnuplot.h.in
	${CMAKE_CURRENT_SOURCE_DIR}/gnuplot.h
)

# The source codes
set(HEADER_FILES
		${CMAKE_CURRENT_SOURCE_DIR}/gnuplot.h	
		Base.h			
		
		HeaderBase.h
		DecoderBase.h 	
		EncoderBase.h
		NetCDFBase.h	

		DWT.h	
		SepDWTHeader.h
		SepDWTPoolBase.h
		SepDWTEncoderPool.h
		SepDWTDecoderPool.h
		WaveletSATEncoder.h
		WaveletSATGPUEncoder.h	

		SATFileNetCDF.h	
		SATFileDecoder.h
		SATFileEncoder.h	
		
		SimpleND.h

		SATSepDWTNetCDF.h	# The class that contains the information for NetCDF format
		SimpleNDFile.h

		SATSepDWTDecoder.h	
		
		utils.h			
	)
	
add_custom_target (
	${PROJECT_NAME}
	SOURCES
		${HEADER_FILES}
)

install(
	FILES ${HEADER_FILES}
	DESTINATION "include")
	

############################################################
# Copyright (c) 2013 Teng-Yok Lee
#
# See the file LICENSE.txt for copying permission.
############################################################
