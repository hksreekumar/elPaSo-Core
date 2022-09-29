#---------------------------------------------------------------------------#
# elPaSo - CMake Project for PETSc
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut f�r Akustik, Technische Universit�t Braunschweig
#---------------------------------------------------------------------------#

# ---- Start point for PETSc libs finding -----------------------------------------------------------#
MESSAGE("> Finding PETSc ...")

IF(${ELPASO_OS_LIN_x86_64} MATCHES FOUND)
    # ---- Set Root Directory -----------------------------------------------------------------------#
    SET(PETSC_ROOT_DIR "${ELPASO_LIB_DIR}/petsc-${PETSC_VERSION}/${PETSC_ARCH}")                # PETSC_DIR
    #SET(PETSC_ROOT_DIR "/home/sreekumar/software/repos/elPaSo_VS/bin")      

    IF(EXISTS ${PETSC_ROOT_DIR})
        # ---- Set each library ---------------------------------------------------------------------#
        IF(LINK_CONAN)
            SET(PETSC_INCLUDE_DIR "${PETSC_ROOT_DIR}/include")  # PETSC_INCLUDE_DIR
        ELSE()
            SET(PETSC_INCLUDE_DIR "${PETSC_ROOT_DIR}/include; ${ELPASO_LIB_DIR}/petsc-${PETSC_VERSION}/include")    # PETSC_INCLUDE_DIR
        ENDIF()
        
        SET(PETSC_LIBRARY_DIR "${PETSC_ROOT_DIR}/lib")          # PETSC_LIBRARY_DIR
        
        # ---- Link and include directories ---------------------------------------------------------#
        LINK_DIRECTORIES(${PETSC_LIBRARY_DIR})
        INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIR})

        # ---- Find each Library --------------------------------------------------------------------#
        SET(HANDLE_PETSC_LIBS)
        IF(LINK_MUMPS)
            FIND_LIBRARY(HANDLE_PETSC_LIB_DMUMPS_${PETSC_ARCH}            libdmumps.a       ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
            FIND_LIBRARY(HANDLE_PETSC_LIB_ZMUMPS_${PETSC_ARCH}            libzmumps.a       ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
            FIND_LIBRARY(HANDLE_PETSC_LIB_CMUMPS_${PETSC_ARCH}            libcmumps.a       ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
            FIND_LIBRARY(HANDLE_PETSC_LIB_SMUMPS_${PETSC_ARCH}            libsmumps.a       ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
            FIND_LIBRARY(HANDLE_PETSC_LIB_MUMPS_COMMON_${PETSC_ARCH}      libmumps_common.a ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
            SET(HANDLE_PETSC_LIBS ${HANDLE_PETSC_LIBS} ${HANDLE_PETSC_LIB_DMUMPS_${PETSC_ARCH}} ${HANDLE_PETSC_LIB_ZMUMPS_${PETSC_ARCH}} ${HANDLE_PETSC_LIB_CMUMPS_${PETSC_ARCH}} ${HANDLE_PETSC_LIB_SMUMPS_${PETSC_ARCH}} ${HANDLE_PETSC_LIB_MUMPS_COMMON_${PETSC_ARCH}})
            ADD_DEFINITIONS(-DPETSC_HAVE_MUMPS)
        ELSE()
            MESSAGE(WARNING "> [EXCLUDED] MUMPS")
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
        
        IF(LINK_CONAN)
            FIND_LIBRARY(HANDLE_PETSC_LIB_PETSC_${PETSC_ARCH}             libpetsc.a          ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
            FIND_LIBRARY(HANDLE_PETSC_LIB_SCALAPACK_${PETSC_ARCH}         libscalapack.a       ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
            FIND_LIBRARY(HANDLE_PETSC_LIB_PARMETIS_${PETSC_ARCH}          libparmetis.a       ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
            FIND_LIBRARY(HANDLE_PETSC_LIB_METIS_${PETSC_ARCH}             libmetis.a          ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
        ELSE()
            FIND_LIBRARY(HANDLE_PETSC_LIB_PETSC_${PETSC_ARCH}             libpetsc.so          ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
            FIND_LIBRARY(HANDLE_PETSC_LIB_SCALAPACK_${PETSC_ARCH}         libscalapack.a       ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
            FIND_LIBRARY(HANDLE_PETSC_LIB_PARMETIS_${PETSC_ARCH}          libparmetis.so       ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
            FIND_LIBRARY(HANDLE_PETSC_LIB_METIS_${PETSC_ARCH}             libmetis.so          ${PETSC_LIBRARY_DIR} NO_DEFAULT_PATH)
        ENDIF()
        
        SET(HANDLE_PETSC_LIBS ${HANDLE_PETSC_LIBS} ${HANDLE_PETSC_LIB_PETSC_${PETSC_ARCH}} ${HANDLE_PETSC_LIB_SCALAPACK_${PETSC_ARCH}} ${HANDLE_PETSC_LIB_PARMETIS_${PETSC_ARCH}} ${HANDLE_PETSC_LIB_METIS_${PETSC_ARCH}})

        # ---- Set FLAGS for PETSc ------------------------------------------------------------------#
        SET(PETSC_COMPILER_FLAGS "")
    ELSE()
        MESSAGE(FATAL_ERROR "> PETSC_ROOT_DIR not found.")
    ENDIF()

ENDIF()

# ---- Display summary ------------------------------------------------------------------------------#
MESSAGE("              PETSC VERSION       : ${PETSC_VERSION}")
MESSAGE("              PETSC ROOT    DIR   : ${PETSC_ROOT_DIR}")
MESSAGE("              PETSC INCLUDE DIR   : ${PETSC_INCLUDE_DIR}")
MESSAGE("              PETSC LIBRARY DIR   : ${PETSC_LIBRARY_DIR}")
MESSAGE("              PETSC LIBS          : ${HANDLE_PETSC_LIBS}")
MESSAGE("              PETSC COMPILER FLAGS: ${PETSC_COMPILER_FLAGS}")

# ---- Assert package handles -----------------------------------------------------------------------#
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSC "> MISSING PETSC LIB HANDLES." HANDLE_PETSC_LIBS)

# ---- Link libraries -------------------------------------------------------------------------------#
LINK_LIBRARIES(${HANDLE_PETSC_LIBS})
