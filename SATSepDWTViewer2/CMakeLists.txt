# Reference:

#The name of the project
PROJECT(SATSepDWTViewer2)

find_package(OpenGL)
find_package(GLUT)

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

set(WITH_FREEGLUT 	FALSE 		CACHE BOOL "Whether FREEGULT is used.")	
set(MYLIB_DIR 		"~/mylib" 	CACHE PATH "The path to MyLib.")	

find_package(ThirdPartyLib)	# ADD-BY-LEETEN 2013/0811

set(WITH_STATIC_LIBS 	TRUE	CACHE BOOL "Whether libraries are statically linked.")
if( WITH_STATIC_LIBS )
	set(WITH_OPENCV		FALSE 	CACHE BOOL "Whether OpenCV is used.")
	if( WITH_OPENCV )
		set(OPENCV_DIR		"C:/OpenCV2.3/vc2008/x64" 	CACHE PATH "The path to OpenCV.")		
		message(WARNING	"For VC projects with static linking: Remember to add the option /NODEFAULTLIB:LIBCMT.lib to CMAKE_EXE_LINKER_FLAGS")

		link_directories (
			"${OPENCV_DIR}/share/OpenCV/3rdparty/lib"
			"${OPENCV_DIR}/lib"
			)
		include_directories (
			"${OPENCV_DIR}/include"
			"${OPENCV_DIR}/include/opencv"
			)
		add_definitions(-DWITH_OPENCV)
	endif( )	

	add_definitions(
		-DWITH_STATIC_LINK 
		-DFREEGLUT_STATIC
		-DGLEW_STATIC
	)
endif()
# ADD-BY-LEETEN 04/04/2012-END
 
if( ${WITH_FREEGLUT} )	
	add_definitions(
		-DUSE_FREEGLUT
	)	
endif()

include_directories(
	${MYLIB_DIR}/include
    # DEL-BY-LEETEN 2013/0811:	${3RDPARTY_DIR}/include	
	${CUDA_TOOLKIT_ROOT_DIR}/include
)

link_directories(
	${MYLIB_DIR}/lib
    # DEL-BY-LEETEN 2013/0811:	${3RDPARTY_DIR}/lib
	${CUDA_TOOLKIT_ROOT_DIR}/lib/Win32
)

# process the shaders
list(APPEND SHADER_NAMES
	ray_integral.frag
	)

foreach(SHADER ${SHADER_NAMES})
	list(APPEND SHADER_FILES	
		${PROJECT_SOURCE_DIR}/${SHADER}
	)
endforeach(SHADER)
	
foreach(SHADER ${SHADER_FILES})
	list(APPEND SHADER_HEADERS
		${SHADER}.h
	)
	file(TO_NATIVE_PATH ${SHADER} SHADER_NATIVE_PATH)
	set(SHADER_OUTPUT_NATIVE_PATH ${SHADER_NATIVE_PATH}.h)
	add_custom_command (
		OUTPUT ${SHADER}.h
		MAIN_DEPENDENCY ${SHADER}
		DEPENDS ${SHADER}
		COMMAND ${MYLIB_DIR}/bin/shader2string ${SHADER_NATIVE_PATH} ${SHADER_OUTPUT_NATIVE_PATH}
	)
endforeach(SHADER)

# The source codes
ADD_EXECUTABLE(${PROJECT_NAME}
	
	SATSepDWTQuery.h

	SATSepDWTBlock.h
	SATSepDWTBlock.cpp

	SATSepDWTBlockTree.h
	
	SATSepDWTViewer_main.cpp
	# The 2D/3D view
	SATSepDWT3DView.h
	SATSepDWT3DView.cpp
	
	# The histogram view
	SATSepDWTHistView.h
	SATSepDWTHistView.cpp
	
	SATSepDWTPCPView.h
	SATSepDWTPCPView.cpp
	
	SATSepDWTQueryView.h
	SATSepDWTQueryView.cpp

	${SHADER_HEADERS}	
	${SHADER_FILES}	
)

# ADD-BY-LEETEN 2013/0811-BEGIN
target_link_libraries(${PROJECT_NAME}	${ThirdPartyLib_LIBRARIES})
# ADD-BY-LEETEN 2013/0811-END

# link OpenGL-related libraries
if(DEFINED WIN32)
	target_link_libraries(${PROJECT_NAME}
		opengl32
		glu32
	    # DEL-BY-LEETEN 2013/0811:	debug	ITKNrrdIOd	optimized	ITKNrrdIO
	)
else()
	target_link_libraries(${PROJECT_NAME}
		GL
		GLU
		GLEW
		glut
	    # DEL-BY-LEETEN 2013/0811-BEGIN
		# glui 
		# netcdf 
		# ITKNrrdIO
	    # DEL-BY-LEETEN 2013/0811-END
	)
endif()

# the needed libraries
if(DEFINED WIN32)	
	target_link_libraries(${PROJECT_NAME}
		debug	libopt_d	optimized	libopt_r
		debug	libgw_d		optimized	libgw_r
		debug	libdvrsuite_d	optimized	libdvrsuite_r
		debug	libdvr2_d	optimized	libdvr2_r
		debug	libclip_d	optimized	libclip_r
		debug	libshader_d	optimized	libshader_r
		debug	libtfw_d	optimized	libtfw_r
	)
	
else()
	target_link_libraries(${PROJECT_NAME}
		libopt
		libgw
		libdvrsuite	
		libdvr2
		libclip
		libshader
		libtfw
		libfbo
		libfps
		# DEL-BY-LEETEN 2013/0811:	lodepng
	)
endif()


############################################################
# Copyright (c) 2013 Teng-Yok Lee
#
# See the file LICENSE.txt for copying permission.
############################################################
