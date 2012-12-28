# Reference:

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#The name of the project
# MOD-BY-LEETEN 10/17/2012-FROM:
	# PROJECT(WaveletSATLib)
# TO:
PROJECT(WaveletSAT)
# MOD-BY-LEETEN 10/17/2012-END

# # DEL-BY-LEETEN 12/28/2012-BEGIN
# set(MYLIB_DIR 		"~/mylib" 	CACHE PATH "The path to MyLib.")	

# include_directories(
	# ${MYLIB_DIR}/include
# )

# link_directories(
	# ${MYLIB_DIR}/lib
# )
# # DEL-BY-LEETEN 12/28/2012-END

# The source codes
set(HEADER_FILES
		HeaderBase.h
		EncoderBase.h
		SATEncoder.h

		DWT.h	
		SepDWTHeader.h
		# DEL-BY-LEETEN 12/25/2012:	SepDWTData.h	# ADD-BY-LEETEN 10/29/2012
		SepDWTPool.h	# ADD-BY-LEETEN 11/11/2012
		WaveletSATEncoder.h
		SATSepDWTEncoder.h
		
		SimpleND.h

		# ADD-BY-LEETEN 12/16/2012-BEGIN
		SATSepDWTNetCDF.h	# The class that contains the information for NetCDF format
		# DEL-BY-LEETEN 12/25/2012:	SATSepDWTFile.h
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
	
