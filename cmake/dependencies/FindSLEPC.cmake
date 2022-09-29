#---------------------------------------------------------------------------#
# elPaSo - CMake Project for SLEPc
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut für Akustik, Technische Universität Braunschweig
#---------------------------------------------------------------------------#

# ---- Start point for SLEPc libs finding -----------------------------------------------------------#
MESSAGE("> Finding SLEPc ...")

IF(${ELPASO_OS_LIN_x86_64} MATCHES FOUND)
    # ---- Set Root Directory -----------------------------------------------------------------------#
    #SET(SLEPC_ROOT_DIR "${ELPASO_LIB_DIR}/slepc-${SLEPC_VERSION}/${ELPASO_COMPILER_OPT}/${SLEPC_ARCH}")                # SLEPC_DIR
    SET(SLEPC_ROOT_DIR "${ELPASO_LIB_DIR}/slepc-${SLEPC_VERSION}/${PETSC_ARCH}")                                        # SLEPC_DIR
    
    IF(EXISTS ${SLEPC_ROOT_DIR})
        # ---- Set each library ---------------------------------------------------------------------#
        SET(SLEPC_INCLUDE_DIR "${SLEPC_ROOT_DIR}/include; ${ELPASO_LIB_DIR}/slepc-${SLEPC_VERSION}/include")    # SLEPC_INCLUDE_DIR
        SET(SLEPC_LIBRARY_DIR "${SLEPC_ROOT_DIR}/lib")        # SLEPC_LIBRARY_DIR
        
        # ---- Link and include directories ---------------------------------------------------------#
        LINK_DIRECTORIES(${SLEPC_LIBRARY_DIR})
        INCLUDE_DIRECTORIES(${SLEPC_INCLUDE_DIR})

        # ---- Find each Library --------------------------------------------------------------------#
        SET(HANDLE_SLEPC_LIBS)

        IF(LINK_CONAN)
            FIND_LIBRARY(HANDLE_SLEPC_LIB_SLPEC_${PETSC_ARCH}         libslepc.a       ${SLEPC_LIBRARY_DIR} NO_DEFAULT_PATH)
        ELSE()
            FIND_LIBRARY(HANDLE_SLEPC_LIB_SLPEC_${PETSC_ARCH}         libslepc.so       ${SLEPC_LIBRARY_DIR} NO_DEFAULT_PATH)
        ENDIF()
        SET(HANDLE_SLEPC_LIBS ${HANDLE_SLEPC_LIBS} ${HANDLE_SLEPC_LIB_SLPEC_${PETSC_ARCH}})

        # ---- Set FLAGS for SLEPC ------------------------------------------------------------------#
        SET(SLEPC_COMPILER_FLAGS "")
    ELSE()
        MESSAGE(FATAL_ERROR "> SLEPC_ROOT_DIR not found.")
    ENDIF()

ENDIF()

# ---- Display summary ------------------------------------------------------------------------------#
MESSAGE("              SLEPC VERSION       : ${SLEPC_VERSION}")
MESSAGE("              SLEPC ROOT    DIR   : ${SLEPC_ROOT_DIR}")
MESSAGE("              SLEPC INCLUDE DIR   : ${SLEPC_INCLUDE_DIR}")
MESSAGE("              SLEPC LIBRARY DIR   : ${SLEPC_LIBRARY_DIR}")
MESSAGE("              SLEPC LIBS          : ${HANDLE_SLEPC_LIBS}")
MESSAGE("              SLEPC COMPILER FLAGS: ${SLEPC_COMPILER_FLAGS}")

# ---- Assert package handles -----------------------------------------------------------------------#
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SLEPC "> MISSING SLEPC LIB HANDLES." HANDLE_SLEPC_LIBS)

# ---- Link libraries -------------------------------------------------------------------------------#
LINK_LIBRARIES(${HANDLE_SLEPC_LIBS})