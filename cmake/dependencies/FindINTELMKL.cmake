#---------------------------------------------------------------------------#
# elPaSo - CMake Project for INTELMKL
#
# 08.11.2021
# Harikrishnan Sreekumar
# Institut f�r Akustik, Technische Universit�t Braunschweig
#---------------------------------------------------------------------------#

# ---- Start point for INTELMKL libs finding -----------------------------------------------------------#
MESSAGE("> Finding INTELMKL ...")

IF(${ELPASO_OS_LIN_x86_64} MATCHES FOUND)
    # ---- Set Root Directory -----------------------------------------------------------------------#
    SET(INTELMKL_ROOT_DIR "${INTEL_DIR}/mkl")                                        # INTELMKL_DIR
    
    IF(EXISTS ${INTELMKL_ROOT_DIR})
        # INTEL LINK ADVISOR
        # LIBS: ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
        # Compiler options:  -DMKL_ILP64  -I"${MKLROOT}/include"

        # ---- Set each library ---------------------------------------------------------------------#
        SET(INTELMKL_INCLUDE_DIR "${INTELMKL_ROOT_DIR};${INTELMKL_ROOT_DIR}/include")    # INTELMKL_INCLUDE_DIR
        SET(INTELMKL_LIBRARY_DIR "${INTELMKL_ROOT_DIR}/lib/intel64/")        # INTELMKL_LIBRARY_DIR
        
        # ---- Link and include directories ---------------------------------------------------------#
        LINK_DIRECTORIES(${INTELMKL_LIBRARY_DIR})
        INCLUDE_DIRECTORIES(${INTELMKL_INCLUDE_DIR})

        # ---- Find each Library --------------------------------------------------------------------#
        SET(HANDLE_INTELMKL_LIBS)
        
        FIND_LIBRARY(HANDLE_INTELMKL_LIB_SCALAPACK          libmkl_scalapack_ilp64.a        ${INTELMKL_LIBRARY_DIR} NO_DEFAULT_PATH)
        FIND_LIBRARY(HANDLE_INTELMKL_LIB_ILP64              libmkl_intel_ilp64.a            ${INTELMKL_LIBRARY_DIR} NO_DEFAULT_PATH)
        FIND_LIBRARY(HANDLE_INTELMKL_LIB_ITHREAD            libmkl_intel_thread.a        ${INTELMKL_LIBRARY_DIR} NO_DEFAULT_PATH)
        FIND_LIBRARY(HANDLE_INTELMKL_LIB_MKLCORE            libmkl_core.so                   ${INTELMKL_LIBRARY_DIR} NO_DEFAULT_PATH)
        FIND_LIBRARY(HANDLE_INTELMKL_LIB_BLACSMPI           libmkl_blacs_intelmpi_ilp64.a                   ${INTELMKL_LIBRARY_DIR} NO_DEFAULT_PATH)
        FIND_LIBRARY(HANDLE_INTELMKL_LIB_LP64               libmkl_intel_lp64.a             ${INTELMKL_LIBRARY_DIR} NO_DEFAULT_PATH)
        #SET(HANDLE_INTELMKL_LIBS ${HANDLE_INTELMKL_LIBS} ${HANDLE_INTELMKL_LIB_ILP64}  ${HANDLE_INTELMKL_LIB_MKLCORE} ${HANDLE_INTELMKL_LIB_ITHREAD}  ${HANDLE_INTELMKL_LIB_BLACSMPI} ${HANDLE_INTELMKL_LIB_LP64})
        SET(HANDLE_INTELMKL_LIBS ${HANDLE_INTELMKL_LIBS} ${HANDLE_INTELMKL_LIB_LP64} ${HANDLE_INTELMKL_LIB_MKLCORE} ${HANDLE_INTELMKL_LIB_ITHREAD})

        # ---- Set FLAGS for INTELMKL ------------------------------------------------------------------#
        #SET(INTELMKL_COMPILER_FLAGS "-liomp5;-lpthread;-lm;-ldl;-DMKL_ILP64")
        SET(INTELMKL_COMPILER_FLAGS "")
    ELSE()
        MESSAGE(FATAL_ERROR "> INTELMKL_ROOT_DIR not found.")
    ENDIF()

ENDIF()

# ---- Display summary ------------------------------------------------------------------------------#
MESSAGE("              INTELMKL VERSION       : ${INTELMKL_VERSION}")
MESSAGE("              INTELMKL ROOT    DIR   : ${INTELMKL_ROOT_DIR}")
MESSAGE("              INTELMKL INCLUDE DIR   : ${INTELMKL_INCLUDE_DIR}")
MESSAGE("              INTELMKL LIBRARY DIR   : ${INTELMKL_LIBRARY_DIR}")
MESSAGE("              INTELMKL LIBS          : ${HANDLE_INTELMKL_LIBS}")
MESSAGE("              INTELMKL COMPILER FLAGS: ${INTELMKL_COMPILER_FLAGS}")

# ---- Assert package handles -----------------------------------------------------------------------#
include(FindPackageHandleStandardArgs)
#find_package_handle_standard_args(INTELMKL "> MISSING INTELMKL LIB HANDLES." HANDLE_INTELMKL_LIBS)

# ---- Link libraries -------------------------------------------------------------------------------#

#ADD_DEFINITIONS(-DMKL_ILP64)
LINK_LIBRARIES("${HANDLE_INTELMKL_LIB_SCALAPACK};${HANDLE_INTELMKL_LIBS}")
LINK_LIBRARIES("-liomp5;-lpthread;-lm;-ldl")
#LINK_LIBRARIES("${INTELMKL_ROOT_DIR}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${INTELMKL_ROOT_DIR}/lib/intel64/libmkl_intel_ilp64.a ${INTELMKL_ROOT_DIR}/lib/intel64/libmkl_intel_thread.a ${INTELMKL_ROOT_DIR}/lib/intel64/libmkl_core.a ${INTELMKL_ROOT_DIR}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group")