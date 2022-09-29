#---------------------------------------------------------------------------#
# elPaSo - CMake Project for INTELMPI
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut f�r Akustik, Technische Universit�t Braunschweig
#---------------------------------------------------------------------------#

# ---- Start point for INTELMPI libs finding -----------------------------------------------------------#
MESSAGE("> Finding INTELMPI ...")

IF(${ELPASO_OS_LIN_x86_64} MATCHES FOUND)
    # ---- Set Root Directory -----------------------------------------------------------------------#
    SET(INTELMPI_ROOT_DIR "${INTEL_DIR}/impi/${INTELMPI_VERSION}/intel64")                                        # INTELMPI_DIR
    
    IF(EXISTS ${INTELMPI_ROOT_DIR})
        # ---- Set each library ---------------------------------------------------------------------#
        SET(INTELMPI_INCLUDE_DIR "${INTELMPI_ROOT_DIR}/include")                                        # INTELMPI_INCLUDE_DIR
        SET(INTELMPI_LIBRARY_DIR "${INTELMPI_ROOT_DIR}/lib" "${INTELMPI_ROOT_DIR}/lib/release")         # INTELMPI_LIBRARY_DIR
        SET(INTELMPI_FABRIC_DIR "${INTELMPI_ROOT_DIR}/libfabric/lib")                                   # INTELMPI_FABRIC_DIR
        SET(INTEL_COMPILER_LIBRARY_DIR "${INTEL_DIR}/compilers_and_libraries/linux/lib/intel64")       # INTEL_COMPILER_LIBRARY_DIR
        
        # ---- Link and include directories ---------------------------------------------------------#
        LINK_DIRECTORIES(${INTELMPI_LIBRARY_DIR} ${INTELMPI_FABRIC_DIR} ${INTEL_COMPILER_LIBRARY_DIR})
        INCLUDE_DIRECTORIES(${INTELMPI_INCLUDE_DIR})

        # ---- Find each Library --------------------------------------------------------------------#
        SET(HANDLE_INTELMPI_LIBS)

        FIND_LIBRARY(HANDLE_INTELMPI_LIB_MPI         libmpi.so       ${INTELMPI_LIBRARY_DIR} NO_DEFAULT_PATH)
        #FIND_LIBRARY(HANDLE_INTELMPI_LIB_FABRIC      libfabric.so    ${INTELMPI_FABRIC_DIR} NO_DEFAULT_PATH)
        SET(HANDLE_INTELMPI_LIBS ${HANDLE_INTELMPI_LIBS} ${HANDLE_INTELMPI_LIB_MPI} ${HANDLE_INTELMPI_LIB_FABRIC})
        
        FIND_LIBRARY(HANDLE_INTELMKL_LIB_FPORT        libifport.so                  ${INTEL_COMPILER_LIBRARY_DIR} NO_DEFAULT_PATH)
        FIND_LIBRARY(HANDLE_INTELMKL_LIB_IMF          libimf.so                     ${INTEL_COMPILER_LIBRARY_DIR} NO_DEFAULT_PATH)
        FIND_LIBRARY(HANDLE_INTELMKL_LIB_INTLC        libintlc.so                   ${INTEL_COMPILER_LIBRARY_DIR} NO_DEFAULT_PATH)
        FIND_LIBRARY(HANDLE_INTELMKL_LIB_VML          libsvml.so                    ${INTEL_COMPILER_LIBRARY_DIR} NO_DEFAULT_PATH)
        FIND_LIBRARY(HANDLE_INTELMKL_LIB_IRC          libirc.so                     ${INTEL_COMPILER_LIBRARY_DIR} NO_DEFAULT_PATH)
        FIND_LIBRARY(HANDLE_INTELMKL_LIB_IRNG         libirng.so                    ${INTEL_COMPILER_LIBRARY_DIR} NO_DEFAULT_PATH)
        SET(HANDLE_INTELMPI_LIBS ${HANDLE_INTELMPI_LIBS} ${HANDLE_INTELMKL_LIB_FPORT} ${HANDLE_INTELMKL_LIB_IMF} ${HANDLE_INTELMKL_LIB_INTLC} ${HANDLE_INTELMKL_LIB_VML} ${HANDLE_INTELMKL_LIB_IRC} ${HANDLE_INTELMKL_LIB_IRNG})
        
        # ---- Set FLAGS for INTELMPI ------------------------------------------------------------------#
        SET(INTELMPI_COMPILER_FLAGS "")
    ELSE()
        MESSAGE(FATAL_ERROR "> INTELMPI_ROOT_DIR not found.")
    ENDIF()

ENDIF()

# ---- Display summary ------------------------------------------------------------------------------#
MESSAGE("              INTELMPI VERSION       : ${INTELMPI_VERSION}")
MESSAGE("              INTELMPI ROOT    DIR   : ${INTELMPI_ROOT_DIR}")
MESSAGE("              INTELMPI INCLUDE DIR   : ${INTELMPI_INCLUDE_DIR}")
MESSAGE("              INTELMPI LIBRARY DIR   : ${INTELMPI_LIBRARY_DIR}")
MESSAGE("              INTELMPI LIBS          : ${HANDLE_INTELMPI_LIBS}")
MESSAGE("              INTELMPI COMPILER FLAGS: ${INTELMPI_COMPILER_FLAGS}")

# ---- Assert package handles -----------------------------------------------------------------------#
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(INTELMPI "> MISSING INTELMPI LIB HANDLES." HANDLE_INTELMPI_LIBS)

# ---- Link libraries -------------------------------------------------------------------------------#
LINK_LIBRARIES(${HANDLE_INTELMPI_LIBS})
