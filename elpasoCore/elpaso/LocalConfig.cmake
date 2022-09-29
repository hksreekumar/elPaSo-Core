#---------------------------------------------------------------------------#
# elPaSo - CMake Project for Module elpaso-LOCAL CONFIG
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
	SET(PETSC_ARCH ${PETSC_ARCH_REAL})
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
	SET(PETSC_ARCH ${PETSC_ARCH_REAL})
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
