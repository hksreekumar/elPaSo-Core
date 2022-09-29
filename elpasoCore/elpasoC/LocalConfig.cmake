#---------------------------------------------------------------------------#
# elPaSo - CMake Project for Module elpasoC-LOCAL CONFIG
#
# 21.09.2021
# Harikrishnan Sreekumar
# Institut für Akustik, Technische Universität Braunschweig
#---------------------------------------------------------------------------#

SET(SPECIFIC_FILES ${ELPASO_SOURCE_DIR}/elpaso/source/main/main.cpp)
SET(ALL_SOURCES ${ALL_SOURCES} ${SPECIFIC_FILES})
SOURCE_GROUP(main FILES ${SPECIFIC_FILES})

# ---- Add relevant .H and .CPP files ----------------------------------#
INCLUDE(${SOURCE_ELPASO}/CMakePackages.cmake)

# ---- Link petsc libraries --------------------------------------------#
IF(LINK_PETSC)
	SET(PETSC_ARCH ${PETSC_ARCH_COMPLEX})
	FIND_PACKAGE(PETSC REQUIRED)
	ADD_DEFINITIONS(-DHAVE_PETSC)
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] PETSC")
ENDIF()	

# ---- Link slepc libraries --------------------------------------------#
IF(LINK_SLEPC)
	FIND_PACKAGE(SLEPC REQUIRED)
	ADD_DEFINITIONS(-DHAVE_SLEPC)
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] SLEPC")
ENDIF()	

# ---- Link arpack libraries -------------------------------------------#
IF(LINK_ARPACK)
	SET(PETSC_ARCH ${PETSC_ARCH_COMPLEX})
	FIND_PACKAGE(ARPACK REQUIRED)
	ADD_DEFINITIONS(-DHAVE_ARPACK)
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] ARPACK")
ENDIF()

# ---- Link mkl libraries ----------------------------------------------#
IF(LINK_MKL)
	IF(ELPASO_COMPILER_INTEL)		
		FIND_PACKAGE(INTELMKL REQUIRED)
		ADD_DEFINITIONS(-DHAVE_INTELMKL)
	ELSEIF(ELPASO_COMPILER_GNU)
		# nothing to add / petsc mkl used
	ENDIF()
	ADD_DEFINITIONS(-DHAVE_MKL)
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] MKL")
ENDIF()

# ---- Link mpi libraries ----------------------------------------------#
IF(LINK_INTELMPI)
	FIND_PACKAGE(INTELMPI REQUIRED)
	ADD_DEFINITIONS(-DHAVE_INTELMPI)
ELSEIF(LINK_OPENMPI)
	FIND_PACKAGE(OPENMPI REQUIRED)
	ADD_DEFINITIONS(-DHAVE_OPENMPI)
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] MPI")
ENDIF()

# ---- Link hdf5 libraries ---------------------------------------------#
IF(LINK_HDF5)
	FIND_PACKAGE(HDF5 REQUIRED)
	ADD_DEFINITIONS(-DHAVE_HDF5)
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] HDF5")
ENDIF()	

# ---- Link xmklio libraries ---------------------------------------------#
IF(LINK_XMLIO)
	FIND_PACKAGE(XMLIO REQUIRED)
	ADD_DEFINITIONS(-DHAVE_XMLIO)
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] XMLIO")
ENDIF()

# ---- Link gzstream libraries ---------------------------------------------#
IF(LINK_GZSTREAM)
	FIND_PACKAGE(GZSTREAM REQUIRED)
	ADD_DEFINITIONS(-DHAVE_GZSTREAM)
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] GZSTREAM")
ENDIF()	

# ---- Clang tidy ------------------------------------------------------#
IF(LINK_CLANGTIDY)		
	set(CMAKE_CXX_CLANG_TIDY
	clang-tidy;
	-header-filter='basics/,elpaso/';
	-checks=-*,clang-analyzer-*,-clang-analyzer-cplusplus*,-clang-analyzer-deadcode.DeadStores,-clang-analyzer-core.UndefinedBinaryOperatorResult,-clang-analyzer-optin.cplusplus.VirtualCall,-clang-analyzer-core.CallAndMessage,-clang-analyzer-core.uninitialized.Branch,-clang-analyzer-security.insecureAPI.strcpy;
	-warnings-as-errors=*;)
	MESSAGE(WARNING "> [INCLUDED] CLANG TIDY")
ELSEIF()
	MESSAGE("> [EXCLUDED] CLANG TIDY")
ENDIF()	

# ---- Set flag for coverage tests ----------------------------------------#
IF(LINK_COV)
	IF(ELPASO_COMPILER_INTEL)
		SET(INTEL_COVERAGE_LINK_FLAGS    "-prof-gen=srcpos")
		SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${INTEL_COVERAGE_LINK_FLAGS}")
		#SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${INTEL_COVERAGE_LINK_FLAGS}")
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