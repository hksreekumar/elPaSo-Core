#---------------------------------------------------------------------------#
# elPaSo - CMake Project
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institute for Acoustics, Technische Universitaet Braunschweig
#
# CMake project rewritten from THE INFAM PROJECT, Marco Schauer, 22.10.2008
#---------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.8)

# ---- Set project name ----------------------------------------------------#
PROJECT("ELPASO")
MESSAGE("> ELPASO base directory: ${ELPASO_SOURCE_DIR}")

# ---- Set source ----------------------------------------------------------#
SET(SOURCE_ELPASO ${ROOT_DIR}/elpasoCore/source)

# ---- Identify target OS --------------------------------------------------#
IF("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
	SET(ELPASOC_CXX_FLAGS_DEBUG   "${ELPASOC_CXX_FLAGS_DEBUG}  " CACHE STRING "ELPASOC_CXX_FLAGS_DEBUG"   FORCE)
	SET(ELPASOC_CXX_FLAGS_RELEASE "${ELPASOC_CXX_FLAGS_RELEASE}" CACHE STRING "ELPASOC_CXX_FLAGS_RELEASE" FORCE)

	SET(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   ${ELPASOC_CXX_FLAGS_DEBUG}")
	SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${ELPASOC_CXX_FLAGS_RELEASE}")

	IF("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "x86_64")
		SET(ELPASO_OS_LIN_x86_64 FOUND)
	ELSE()
		SET(ELPASO_OS_LIN NOTFOUND)
		MESSAGE(FATAL_ERROR "> Unrecognized OS ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
	ENDIF()
ELSE()
	SET(ELPASO_OS NOTFOUND)
	MESSAGE(FATAL_ERROR "> Unrecognized OS ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
ENDIF()
MESSAGE("> Build OS set to: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")

# ---- Set default compiler ------------------------------------------------#
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)

# ---- Identify target Compiler --------------------------------------------#
MESSAGE("> Identified compilers")
MESSAGE("              C   compiler | ID: ${CMAKE_C_COMPILER_ID} | VERSION: ${CMAKE_C_COMPILER_VERSION} | Path: ${CMAKE_C_COMPILER}")
MESSAGE("              CXX compiler | ID: ${CMAKE_CXX_COMPILER_ID} | VERSION: ${CMAKE_CXX_COMPILER_VERSION} | Path: ${CMAKE_CXX_COMPILER}")

IF("${CMAKE_C_COMPILER_ID}" MATCHES "Intel" AND "${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
	SET(ELPASO_COMPILER_INTEL FOUND)
	STRING(TOLOWER ${CMAKE_CXX_COMPILER_ID} ELPASO_COMPILER_ID)
ELSEIF("${CMAKE_C_COMPILER_ID}" MATCHES "GNU" AND "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
	SET(ELPASO_COMPILER_GNU FOUND)
	STRING(TOLOWER ${CMAKE_CXX_COMPILER_ID} ELPASO_COMPILER_ID)
ELSE()
	SET(ELPASO_COMPILER NOTFOUND)
	MESSAGE(FATAL_ERROR "> Unrecognized compiler")
ENDIF()
MESSAGE("> Compiler set to: ${ELPASO_COMPILER_ID}")

# ---- Compiler Option ----------------------------------------------------#
SET(ELPASO_COMPILER_OPT "${ELPASO_COMPILER_ID}-opt")

# ---- Set basic compiler flags -------------------------------------------#
SET(CMAKE_C_FLAGS "")
SET(CMAKE_CXX_FLAGS "")

# ---- Pick dependencies --------------------------------------------------#
INCLUDE(${ELPASO_SOURCE_DIR}/cmake/macro/CMakeElpasoDependencies.cmake)
SET_ELPASO_DEPENDENCIES()

# ---- Add 3rdParty dependencies
INCLUDE(CMakeElpaso3rdPartyDependencies.cmake)
SET_ELPASO_3RDPARTY_DEPENDENCIES()

# ---- Set dependency module path ------------------------------------------#
SET(CMAKE_MODULE_PATH "${ELPASO_SOURCE_DIR}/cmake/dependencies")

# ---- Set MPI according to which compiler ---------------------------------#
IF(ELPASO_COMPILER_INTEL)
	OPTION(LINK_INTELMPI	"MPI from Intel"					ON)
ELSEIF(ELPASO_COMPILER_GNU)
	OPTION(LINK_OPENMPI		"MPI from OpenMPI"					ON)
ENDIF()

# ---- Set flag for coverage tests ----------------------------------------#
IF(LINK_COV)
	IF(ELPASO_COMPILER_INTEL)
		SET(INTEL_COVERAGE_LINK_FLAGS    "-prof-gen=srcpos")
		SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${INTEL_COVERAGE_LINK_FLAGS}")
		SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${INTEL_COVERAGE_LINK_FLAGS}")
		MESSAGE("> Intel codecov activated")
	ELSEIF(ELPASO_COMPILER_GNU)
		SET(GCC_COVERAGE_COMPILE_FLAGS "-g -fprofile-arcs -ftest-coverage")
		SET(GCC_COVERAGE_LINK_FLAGS    "-lgcov")
		SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
		SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
		MESSAGE("> GLCOV activated")
	ENDIF()
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] COV")
ENDIF()

# ---- CMake individual projects -------------------------------------------#
IF(LINK_ELPASO)
	ADD_SUBDIRECTORY(./elpasoCore/elpaso)
ENDIF()

IF(LINK_ELPASOC)
	ADD_SUBDIRECTORY(./elpasoCore/elpasoC)
ENDIF()

IF(LINK_ELPASOT)
	ADD_SUBDIRECTORY(./elpasoCore/elpasoT)
ENDIF()

