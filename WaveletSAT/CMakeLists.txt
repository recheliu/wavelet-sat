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
		EncoderBase.h
		SATEncoder.h

		DWT.h	
		SepDWTHeader.h
		SepDWTPool.h	# ADD-BY-LEETEN 11/11/2012
		WaveletSATEncoder.h
		SATSepDWTEncoder.h
		
		SimpleND.h

		# ADD-BY-LEETEN 12/16/2012-BEGIN
		SATSepDWTNetCDF.h	# The class that contains the information for NetCDF format
		SimpleNDFile.h
		# ADD-BY-LEETEN 12/16/2012-END

		SATSepDWTOutOfCore.h	# ADD-BY-LEETEN 12/25/2012
		
		utils.h			# ADD-BY-LEETEN 10/29/2012
	)
	
# ADD-BY-LEETEN 10/28/2012-BEGIN
if( WITH_NETCDF )
	list(APPEND HEADER_FILES
		SATFileEncoder.h	
	)
endif()
# ADD-BY-LEETEN 10/28/2012-END
			
add_custom_target (
	${PROJECT_NAME}
	SOURCES
		${HEADER_FILES}
)

install(
	FILES ${HEADER_FILES}
	DESTINATION "include")
	
