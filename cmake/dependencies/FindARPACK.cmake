#---------------------------------------------------------------------------#
# elPaSo - CMake Project for ARPACK
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut f�r Akustik, Technische Universit�t Braunschweig
#---------------------------------------------------------------------------#

# ---- Start point for ARPACK libs finding -----------------------------------------------------------#
MESSAGE("> Finding ARPACK ...")

IF(${ELPASO_OS_LIN_x86_64} MATCHES FOUND)
    # ---- Set Root Directory -----------------------------------------------------------------------#
    SET(ARPACK_ROOT_DIR "${ELPASO_LIB_DIR}/ARPACK")             # ARPACK_DIR
    
    IF(EXISTS ${ARPACK_ROOT_DIR})
        # ---- Set each library ---------------------------------------------------------------------#
        SET(ARPACK_INCLUDE_DIR "${ARPACK_ROOT_DIR}")    # ARPACK_INCLUDE_DIR
        SET(ARPACK_LIBRARY_DIR "${ARPACK_ROOT_DIR}/${PETSC_ARCH}/lib")        # ARPACK_LIBRARY_DIR
        
        # ---- Link and include directories ---------------------------------------------------------#
        LINK_DIRECTORIES(${ARPACK_LIBRARY_DIR})
        INCLUDE_DIRECTORIES(${ARPACK_INCLUDE_DIR})

        # ---- Find each Library --------------------------------------------------------------------#
        SET(HANDLE_ARPACK_LIBS)

        FIND_LIBRARY(HANDLE_ARPACK_LIB_ARPACK_${PETSC_ARCH}         libarpack_SUN4.a            ${ARPACK_LIBRARY_DIR} NO_DEFAULT_PATH)
        FIND_LIBRARY(HANDLE_ARPACK_LIB_PARPACK_${PETSC_ARCH}        libparpack_MPI-SUN4.a       ${ARPACK_LIBRARY_DIR} NO_DEFAULT_PATH)
        SET(HANDLE_ARPACK_LIBS ${HANDLE_ARPACK_LIBS} ${HANDLE_ARPACK_LIB_ARPACK_${PETSC_ARCH}} ${HANDLE_ARPACK_LIB_PARPACK_${PETSC_ARCH}})

        # ---- Set FLAGS for ARPACK ------------------------------------------------------------------#
        SET(ARPACK_COMPILER_FLAGS "")
    ELSE()
        MESSAGE(FATAL_ERROR "> ARPACK_ROOT_DIR not found.")
    ENDIF()

ENDIF()

# ---- Display summary ------------------------------------------------------------------------------#
MESSAGE("              ARPACK ROOT    DIR   : ${ARPACK_ROOT_DIR}")
MESSAGE("              ARPACK INCLUDE DIR   : ${ARPACK_INCLUDE_DIR}")
MESSAGE("              ARPACK LIBRARY DIR   : ${ARPACK_LIBRARY_DIR}")
MESSAGE("              ARPACK LIBS          : ${HANDLE_ARPACK_LIBS}")
MESSAGE("              ARPACK COMPILER FLAGS: ${ARPACK_COMPILER_FLAGS}")

# ---- Assert package handles -----------------------------------------------------------------------#
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARPACK "> MISSING ARPACK LIB HANDLES." HANDLE_ARPACK_LIBS)

# ---- Link libraries -------------------------------------------------------------------------------#
LINK_LIBRARIES(${HANDLE_ARPACK_LIBS})