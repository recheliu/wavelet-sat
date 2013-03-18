# Reference:

#The name of the project
PROJECT(SATSepDWTViewer)

find_package(OpenGL)
find_package(GLUT)

#The CMake Minimum version that is required. The FindCUDA script
#is distributed since version 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

set(WITH_FREEGLUT 	FALSE 		CACHE BOOL "Whether FREEGULT is used.")	
set(MYLIB_DIR 		"~/mylib" 	CACHE PATH "The path to MyLib.")	

set(WITH_STATIC_LIBS 	TRUE	CACHE BOOL "Whether libraries are statically linked.")
if( WITH_STATIC_LIBS )
	set(WITH_OPENCV		FALSE 	CACHE BOOL "Whether OpenCV is used.")
	if( WITH_OPENCV )
		set(OPENCV_DIR		"C:/OpenCV2.3/vc2008/x64" 	CACHE PATH "The path to OpenCV.")		
		message(WARNING	"For VC projects with static linking: Remember to add the option /NODEFAULTLIB:LIBCMT.lib to CMAKE_EXE_LINKER_FLAGS")

		link_directories (
			"${OPENCV_DIR}/share/OpenCV/3rdparty/lib"
			"${OPENCV_DIR}/lib"	# ADD-BY-LEETEN 01/19/2013
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
	${3RDPARTY_DIR}/include	
	${CUDA_TOOLKIT_ROOT_DIR}/include
)

link_directories(
	${MYLIB_DIR}/lib
	${3RDPARTY_DIR}/lib
	${CUDA_TOOLKIT_ROOT_DIR}/lib/Win32
)

# process the shaders
list(APPEND SHADER_NAMES
	ray_integral.frag	# ADD-BY-LEETEN 02/05/2013
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
	# The header that define the 
	SATSepDWT.h	# 
	
	SATSepDWTViewer_main.cpp
	# The 2D/3D view
	SATSepDWT3DView.h
	SATSepDWT3DView.cpp
	
	# The histogram view
	SATSepDWTHistView.h
	SATSepDWTHistView.cpp
	
	# ADD-BY-LEETEN 02/14/2013-BEGIN
	SATSepDWTPCPView.h
	SATSepDWTPCPView.cpp
	# ADD-BY-LEETEN 02/14/2013-END
	
	# ADD-BY-LEETEN 03/17/2013-BEGIN
	SATSepDWTQueryView.h
	SATSepDWTQueryView.cpp
	# ADD-BY-LEETEN 03/17/2013-END
	
	${SHADER_HEADERS}	
	${SHADER_FILES}	
)

# link OpenGL-related libraries
if(DEFINED WIN32)
	target_link_libraries(${PROJECT_NAME}
		opengl32
		glu32
		debug	ITKNrrdIOd	optimized	ITKNrrdIO
	)
else()
	target_link_libraries(${PROJECT_NAME}
		GL
		GLU
		GLEW
		glut
		ITKNrrdIO
	)
endif()

# the needed libraries
if(DEFINED WIN32)	
	target_link_libraries(${PROJECT_NAME}
		debug	libopt_d	optimized	libopt_r
		debug	libgw_d		optimized	libgw_r
		debug	libdvrsuite_d	optimized	libdvrsuite_r	# ADD-BY-LEETEN 02/10/2013
		debug	libdvr2_d	optimized	libdvr2_r
		debug	libclip_d	optimized	libclip_r
		debug	libshader_d	optimized	libshader_r
		debug	libtfw_d	optimized	libtfw_r	# ADD-BY-LEETEN 02/05/2013
	)
	
else()
	target_link_libraries(${PROJECT_NAME}
		libopt
		libgw
		libdvrsuite	# ADD-BY-LEETEN 02/10/2013
		libdvr2
		libclip
		libshader
		libtfw	# ADD-BY-LEETEN 02/05/2013
	)
endif()

