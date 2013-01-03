# Reference:

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
PROJECT(WaveletSAT)

# The source codes
set(HEADER_FILES
		Base.h			# ADD-BY-LEETEN 12/30/2012
		
		HeaderBase.h
		DecoderBase.h 	# ADD-BY-LEETEN 01/02/2013
		EncoderBase.h
		NetCDFBase.h	# ADD-BY-LEETEN 01/02/2013
		# DEL-BY-LEETEN 01/02/2013:	SATEncoder.h

		DWT.h	
		SepDWTHeader.h
		SepDWTPool.h	# ADD-BY-LEETEN 11/11/2012
		WaveletSATEncoder.h
		# DEL-BY-LEETEN 01/02/2013:	SATSepDWTEncoder.h

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

		# MOD-BY-LEETEN 01/02/2013-FROM:
			# SATSepDWTOutOfCore.h	# ADD-BY-LEETEN 12/25/2012
		# TO:
		SATSepDWTDecoder.h	# ADD-BY-LEETEN 12/25/2012
		# MOD-BY-LEETEN 01/02/2013-END
		
		utils.h			# ADD-BY-LEETEN 10/29/2012
	)
	
# # DEL-BY-LEETEN 01/02/2013-BEGIN
# # ADD-BY-LEETEN 10/28/2012-BEGIN
# if( WITH_NETCDF )
	# list(APPEND HEADER_FILES
		# # ADD-BY-LEETEN 12/16/2012-BEGIN
		# SATFileNetCDF.h	
		# SATFileDecoder.h
		# # ADD-BY-LEETEN 12/16/2012-END
		# SATFileEncoder.h	
	# )
# endif()
# # ADD-BY-LEETEN 10/28/2012-END
# # DEL-BY-LEETEN 01/02/2013-END
			
add_custom_target (
	${PROJECT_NAME}
	SOURCES
		${HEADER_FILES}
)

install(
	FILES ${HEADER_FILES}
	DESTINATION "include")
	
