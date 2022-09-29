#---------------------------------------------------------------------------#
# elPaSo - CMake Project for Module elpasoC
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut für Akustik, Technische Universität Braunschweig
#
# CMake project rewritten from THE INFAM PROJECT, Marco Schauer, 22.10.2008
#---------------------------------------------------------------------------#

MESSAGE("> Building elPaSo complex module: elpasoC")

cmake_minimum_required(VERSION 3.8)
PROJECT(elpasoC)

IF(${ELPASO_OS_LIN_x86_64} MATCHES FOUND)
	MESSAGE("> Build initiated for Linux")
	# ---- Basic -----------------------------------------------------------#
	SET(EXECUTABLE_NAME elpasoC)
	SET(EXECUTABLE_OUTPUT_PATH ${ELPASO_SOURCE_DIR}/bin)
	SET(CMAKE_BUILD_TYPE Release)

	SET(ROOT_DIR ${ELPASO_SOURCE_DIR})
	SET(SOURCE_ELPASO ${ROOT_DIR}/elpaso/source)
	
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

	INCLUDE(${ELPASO_SOURCE_DIR}/cmake/macro/CMakeElpasoSetGitRepository.cmake)
	SET_ELPASO_GIT()
	
	ADD_DEFINITIONS(-D_MT -DNO_BUILDINFO_H -D_CONSOLE -D_USE_MATH_DEFINES -DPETSC_USE_COMPLEX -DPETSC_CLANGUAGE_CXX)
	SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${GENERATED_FILES}")
ELSE()
	MESSAGE(FATAL_ERROR "> BUILD CANNOT PROCEED. CHECK OS.")
ENDIF()

# Set output path for the new conan routine
SET(CMAKE_CURRENT_BINARY_DIR ${EXECUTABLE_OUTPUT_PATH})

ADD_EXECUTABLE( ${EXECUTABLE_NAME} ${ALL_SOURCES} )
target_link_libraries(${EXECUTABLE_NAME} ${CONAN_TARGETS})