#---------------------------------------------------------------------------#
# elPaSo - CMake Project for Module elpasoT-LOCAL CONFIG
#
# 21.09.2021
# Harikrishnan Sreekumar
# Institut f�r Akustik, Technische Universit�t Braunschweig
#---------------------------------------------------------------------------#

IF(LINK_CLANGTIDY)
	set(CMAKE_CXX_CLANG_TIDY
	clang-tidy;
	-header-filter=.;
	-checks=*;)
	MESSAGE(WARNING "> [INCLUDED] CLANG TIDY")
ELSEIF()
	MESSAGE("> [EXCLUDED] CLANG TIDY")
ENDIF()

SET(SPECIFIC_FILES ${ELPASO_SOURCE_DIR}/elpasoCore/test/testmain.cpp)
SET(ALL_SOURCES ${ALL_SOURCES} ${SPECIFIC_FILES})

SOURCE_GROUP(main FILES ${SPECIFIC_FILES})

INCLUDE(${SOURCE_ELPASO}/CMakePackages.cmake)

# ---- Link gtest libraries
IF(LINK_GTEST)
	FIND_PACKAGE(GTEST REQUIRED)
	ADD_DEFINITIONS(-DHAVE_GTEST)
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] GTEST")
ENDIF()

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