#---------------------------------------------------------------------------#
# elPaSo - CMake Project for Module elpasoT-CONAN CONFIG
#
# 21.09.2021
# Harikrishnan Sreekumar
# Institut f�r Akustik, Technische Universit�t Braunschweig
#---------------------------------------------------------------------------#

# Set the build directory acc. to conan
SET(CMAKE_PROJECT_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})

# For the CONAN generator cmake_find_package, module path is changed
list(APPEND CMAKE_MODULE_PATH ${CMAKE_PROJECT_BUILD_DIR})
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_PROJECT_BUILD_DIR})

# Configure required packages and install
IF(ELPASO_COMPILER_GNU)
	MESSAGE("> Conan installation for gnu compilers...")
	conan_cmake_configure(REQUIRES petsc-complex/${PETSC_VERSION}@ina+elpaso/stable 
				REQUIRES slepc-complex/${SLEPC_VERSION}@ina+elpaso/stable
				REQUIRES arpack-complex/2.1@ina+elpaso/stable
				REQUIRES hdf5/${HDF5_VERSION}@ina+elpaso/stable
				REQUIRES openmpi/${OPENMPI_VERSION}@ina+elpaso/stable
				REQUIRES 
				GENERATORS cmake_find_package
				BASIC_SETUP SETTINGS compiler=gcc SETTINGS compiler.libcxx=libstdc++11 SETTINGS compiler.version=8.4
                IMPORTS "include, *.h* -> ${EXECUTABLE_OUTPUT_PATH}/include"
				IMPORTS "lib, *.so* -> ${EXECUTABLE_OUTPUT_PATH}/lib"
				IMPORTS "lib, *.a* -> ${EXECUTABLE_OUTPUT_PATH}/lib"
				IMPORTS "bin, * -> ${EXECUTABLE_OUTPUT_PATH}/bin"
				IMPORTS "share, * -> ${EXECUTABLE_OUTPUT_PATH}/share")
ELSEIF(ELPASO_COMPILER_INTEL)
	MESSAGE("> Conan installation for intel compilers...")
	conan_cmake_configure(REQUIRES petsc-complex/${PETSC_VERSION}@ina+elpaso/stable 
				REQUIRES slepc-complex/${SLEPC_VERSION}@ina+elpaso/stable
				REQUIRES arpack-complex/2.1@ina+elpaso/stable
				REQUIRES hdf5/${HDF5_VERSION}@ina+elpaso/stable
				REQUIRES 
				GENERATORS cmake_find_package
				BASIC_SETUP SETTINGS compiler=intel SETTINGS compiler.libcxx=libstdc++11 SETTINGS compiler.version=19.1)
ELSE()
	MESSAGE(FATAL_ERROR "> No compiler setting for conan installation")
ENDIF()

conan_cmake_autodetect(settings)
conan_cmake_install(PATH_OR_REFERENCE .
                BUILD missing
                SETTINGS ${settings})

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
	FIND_PACKAGE(petsc-complex REQUIRED)
	list(APPEND CONAN_TARGETS petsc-complex::petsc-complex)
	ADD_DEFINITIONS(-DHAVE_PETSC)
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] PETSC")
ENDIF()
	
# ---- Link slepc libraries --------------------------------------------#
IF(LINK_SLEPC)
	FIND_PACKAGE(slepc-complex REQUIRED)
	list(APPEND CONAN_TARGETS slepc-complex::slepc-complex)
	ADD_DEFINITIONS(-DHAVE_SLEPC)
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] SLEPC")
ENDIF()
	
# ---- Link arpack libraries -------------------------------------------#
IF(LINK_ARPACK)
	SET(PETSC_ARCH ${PETSC_ARCH_COMPLEX})
	FIND_PACKAGE(arpack-complex REQUIRED)
	list(APPEND CONAN_TARGETS arpack-complex::arpack-complex)
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
	FIND_PACKAGE(openmpi REQUIRED)
	list(APPEND CONAN_TARGETS openmpi::openmpi)
	ADD_DEFINITIONS(-DHAVE_OPENMPI)
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] MPI")
ENDIF()

# ---- Link hdf5 libraries ---------------------------------------------#
IF(LINK_HDF5)
	FIND_PACKAGE(hdf5 REQUIRED)
	list(APPEND CONAN_TARGETS hdf5::hdf5)
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
	
# ---- Set flag for coverage tests ----------------------------------------#
IF(LINK_COV)
	IF(ELPASO_COMPILER_INTEL)
		SET(INTEL_COVERAGE_LINK_FLAGS    "-prof-gen=srcpos")
		SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${INTEL_COVERAGE_LINK_FLAGS}")
		SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${INTEL_COVERAGE_LINK_FLAGS}")
		MESSAGE("> Intel codecov activated")
	ELSEIF(ELPASO_COMPILER_GNU)
		SET(GCC_COVERAGE_COMPILE_FLAGS "-g -fprofile-arcs -ftest-coverage")
		SET(GCC_COVERAGE_LINK_FLAGS    "-g -fprofile-arcs -ftest-coverage")
		SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
		SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
		MESSAGE("> GLCOV activated")
	ENDIF()
ELSE()
	MESSAGE(WARNING "> [EXCLUDED] COV")
ENDIF()

IF(ELPASO_COMPILER_INTEL)
    IF(LINK_PARDISO)
        # nothing to include
        SET(HANDLE_PETSC_LIBS ${HANDLE_PETSC_LIBS})
        ADD_DEFINITIONS(-DPETSC_HAVE_PARDISO)
    ELSE()
        MESSAGE(WARNING "> [EXCLUDED] PARDISO")
    ENDIF()
ENDIF()