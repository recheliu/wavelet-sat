# Reference:

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(WaveletSAT)

# ADD-BY-LEETEN 01/23/2013-BEGIN
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
# ADD-BY-LEETEN 01/23/2013-END

# The source codes
set(HEADER_FILES
		${CMAKE_CURRENT_SOURCE_DIR}/gnuplot.h	# ADD-BY-LEETEN 01/23/2013
		Base.h			# ADD-BY-LEETEN 12/30/2012
		
		HeaderBase.h
		DecoderBase.h 	# ADD-BY-LEETEN 01/02/2013
		EncoderBase.h
		NetCDFBase.h	# ADD-BY-LEETEN 01/02/2013

		DWT.h	
		SepDWTHeader.h
		# MOD-BY-LEETEN 2013/10/27:	SepDWTPool.h	# ADD-BY-LEETEN 11/11/2012
		SepDWTPoolBase.h
		SepDWTEncoderPool.h
		SepDWTDecoderPool.h
		# MOD-BY-LEETEN 2013/10/27-END
		WaveletSATEncoder.h
		WaveletSATGPUEncoder.h	# ADD-BY-LEETEN 01/10/2012

		# ADD-BY-LEETEN 01/02/2013-BEGIN
		SATFileNetCDF.h	
		SATFileDecoder.h
		SATFileEncoder.h	
		# ADD-BY-LEETEN 01/02/2013-END
		
		SimpleND.h

		# ADD-BY-LEETEN 12/16/2012-BEGIN
		SATSepDWTNetCDF.h	# The class that contains the information for NetCDF format
		SimpleNDFile.h
		# ADD-BY-LEETEN 12/16/2012-END

		SATSepDWTDecoder.h	# ADD-BY-LEETEN 12/25/2012
		
		utils.h			# ADD-BY-LEETEN 10/29/2012
	)
	
add_custom_target (
	${PROJECT_NAME}
	SOURCES
		${HEADER_FILES}
)

install(
	FILES ${HEADER_FILES}
	DESTINATION "include")
	
