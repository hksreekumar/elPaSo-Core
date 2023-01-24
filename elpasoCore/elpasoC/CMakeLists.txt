#---------------------------------------------------------------------------#
# elPaSo - CMake Project for Module elpasoC
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut f�r Akustik, Technische Universit�t Braunschweig
#
# CMake project rewritten from THE INFAM PROJECT, Marco Schauer, 22.10.2008
#---------------------------------------------------------------------------#

MESSAGE("> Building elPaSo complex module: elpasoC")

SET(ELPASO_EXEC_NAME "elpasoC")
SET(ELPASO_DLIB_NAME "elpasoCore-cxx-complex-o")

cmake_minimum_required(VERSION 3.8)
IF(NOT GEN_DLIB)
	SET(ELPASO_PROJECT_NAME ${ELPASO_EXEC_NAME})	
ELSE()
	SET(ELPASO_PROJECT_NAME ${ELPASO_DLIB_NAME})
ENDIF()
PROJECT(${ELPASO_PROJECT_NAME})

IF(${ELPASO_OS_LIN_x86_64} MATCHES FOUND)
	MESSAGE("> Build initiated for Linux")
	# ---- Basic -----------------------------------------------------------#
	SET(EXECUTABLE_NAME elpasoC)
	SET(EXECUTABLE_OUTPUT_PATH ${ELPASO_SOURCE_DIR}/bin)
	SET(CMAKE_BUILD_TYPE Release)

	SET(ROOT_DIR ${ELPASO_SOURCE_DIR})
	SET(SOURCE_ELPASO ${ROOT_DIR}/elpasoCore/source)
	
	IF(LINK_CONAN)
		MESSAGE(STATUS "ELPASO INSTALLATION USING CONAN REPOS")
		# Download automatically, you can also just copy the conan.cmake file
		if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
		   message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
		   file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/master/conan.cmake"
						  "${CMAKE_BINARY_DIR}/conan.cmake")
		endif()
		include(${CMAKE_BINARY_DIR}/conan.cmake)

		include(./ConanConfig.cmake)
	ELSE()
		MESSAGE(STATUS "ELPASO INSTALLATION USING LOCAL REPOS")
		include(./LocalConfig.cmake)
	ENDIF()

	# ---- Clang tidy ------------------------------------------------------#
	IF(LINK_CLANGTIDY)		
		set(CMAKE_CXX_CLANG_TIDY
		clang-tidy;
		-header-filter='elpasoCore/';
		-checks=-*,clang-analyzer-*,-clang-analyzer-cplusplus*,-clang-analyzer-deadcode.DeadStores,-clang-analyzer-core.UndefinedBinaryOperatorResult,-clang-analyzer-optin.cplusplus.VirtualCall,-clang-analyzer-core.CallAndMessage,-clang-analyzer-core.uninitialized.Branch,-clang-analyzer-security.insecureAPI.strcpy;
		-warnings-as-errors=*;)
		MESSAGE(WARNING "> [INCLUDED] CLANG TIDY")
	ELSEIF()
		MESSAGE("> [EXCLUDED] CLANG TIDY")
	ENDIF()	

	INCLUDE(${ELPASO_SOURCE_DIR}/cmake/macro/CMakeElpasoSetGitRepository.cmake)
	SET_ELPASO_GIT()
	
	ADD_DEFINITIONS(-D_MT -DNO_BUILDINFO_H -D_CONSOLE -D_USE_MATH_DEFINES -DPETSC_USE_COMPLEX -DPETSC_CLANGUAGE_CXX)
	SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${GENERATED_FILES}")
ELSE()
	MESSAGE(FATAL_ERROR "> BUILD CANNOT PROCEED. CHECK OS.")
ENDIF()

# Set output path for the new conan routine
SET(CMAKE_CURRENT_BINARY_DIR ${EXECUTABLE_OUTPUT_PATH})

#MESSAGE( "> ALL SOURCES ${ALL_SOURCES}")
# Generating executable
IF(NOT GEN_DLIB)
	MESSAGE( "> EXECUTABLE GENERATION ACTIVATED")
	ADD_EXECUTABLE( ${EXECUTABLE_NAME} ${ALL_SOURCES} )
	TARGET_LINK_LIBRARIES(${EXECUTABLE_NAME} ${CONAN_TARGETS})
ENDIF()

# Generating library
IF(GEN_DLIB)
	SET(CMAKE_INSTALL_PREFIX 	${ELPASO_SOURCE_DIR}/elpasoCore-${ELPASO_COMPILER_ID} CACHE PATH "LIB INSTALLATION" FORCE)
	SET(DEFAULT_BUILD_TYPE 		"Release")
	SET(ALL_HEADERS 			${ALL_SOURCES})
	SET(ALL_CPP 				${ALL_SOURCES})
	
	INCLUDE(GNUInstallDirs)

	LIST(FILTER ALL_HEADERS INCLUDE REGEX ".h") 
	LIST(FILTER ALL_CPP 	INCLUDE REGEX ".cpp") 
	
	ADD_LIBRARY(${PROJECT_NAME} SHARED ${ALL_CPP})
	#TARGET_INCLUDE_DIRECTORIES(${PROJECT_NAME} PRIVATE ${SOURCE_ELPASO})
	#MESSAGE( "> ALL SOURCES ${ALL_SOURCES}")
	#MESSAGE( "> ALL SOURCES ${ALL_CPP}")
	#SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${ALL_HEADERS}")

	INSTALL(TARGETS ${PROJECT_NAME}
			EXPORT "${PROJECT_NAME}Config"
			LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
			PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/elpasoCore
			)
	INSTALL(EXPORT "${PROJECT_NAME}Config"
			FILE "${PROJECT_NAME}Config.cmake"
			DESTINATION share/cmake
			)

	FILE(GLOB SOURCE_FOLDER_LIST LIST_DIRECTORIES true ${SOURCE_ELPASO}/*)
	FOREACH(SOURCE_FOLDER ${SOURCE_FOLDER_LIST})
		IF(IS_DIRECTORY ${SOURCE_FOLDER})
			INSTALL(DIRECTORY ${SOURCE_FOLDER}
					DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/elpasoCore
					FILES_MATCHING
					PATTERN "*.h"
				)
		ENDIF()
	ENDFOREACH()

	TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${CONAN_TARGETS})
	MESSAGE( "> DYNAMIC LIBRARY GENERATION ACTIVATED")
	MESSAGE( "              TARGET NAME : ${PROJECT_NAME}")
	MESSAGE( "              INSTALL PATH: ${CMAKE_INSTALL_PREFIX}")
ENDIF()