#---------------------------------------------------------------------------#
# elPaSo - CMake Project for OPENMPI
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut für Akustik, Technische Universität Braunschweig
#---------------------------------------------------------------------------#

# ---- Start point for OPENMPI libs finding -----------------------------------------------------------#
MESSAGE("> Finding OPENMPI ...")

IF(${ELPASO_OS_LIN_x86_64} MATCHES FOUND)
    # ---- Set Root Directory -----------------------------------------------------------------------#
    SET(OPENMPI_ROOT_DIR "${ELPASO_LIB_DIR}/openmpi-${OPENMPI_VERSION}/${ELPASO_COMPILER_OPT}")                                        # OPENMPI_DIR
    
    IF(EXISTS ${OPENMPI_ROOT_DIR})
        # ---- Set each library ---------------------------------------------------------------------#
        SET(OPENMPI_INCLUDE_DIR "${OPENMPI_ROOT_DIR}/include")    # OPENMPI_INCLUDE_DIR
        SET(OPENMPI_LIBRARY_DIR "${OPENMPI_ROOT_DIR}/lib")        # OPENMPI_LIBRARY_DIR
        
        # ---- Link and include directories ---------------------------------------------------------#
        LINK_DIRECTORIES(${OPENMPI_LIBRARY_DIR})
        INCLUDE_DIRECTORIES(${OPENMPI_INCLUDE_DIR})

        # ---- Find each Library --------------------------------------------------------------------#
        SET(HANDLE_OPENMPI_LIBS)
        
	    #IF(LINK_CONAN)
        #    FIND_LIBRARY(HANDLE_OPENMPI_LIB_MPI         libmpi.a       ${OPENMPI_LIBRARY_DIR} NO_DEFAULT_PATH)
        #ELSE()
            FIND_LIBRARY(HANDLE_OPENMPI_LIB_MPI         libmpi.so      ${OPENMPI_LIBRARY_DIR} NO_DEFAULT_PATH)
        #ENDIF()
        SET(HANDLE_OPENMPI_LIBS ${HANDLE_OPENMPI_LIBS} ${HANDLE_OPENMPI_LIB_MPI})

        # ---- Set FLAGS for OPENMPI ------------------------------------------------------------------#
        SET(OPENMPI_COMPILER_FLAGS "")
    ELSE()
        MESSAGE(FATAL_ERROR "> OPENMPI_ROOT_DIR not found.")
    ENDIF()

ENDIF()

# ---- Display summary ------------------------------------------------------------------------------#
MESSAGE("              OPENMPI VERSION       : ${OPENMPI_VERSION}")
MESSAGE("              OPENMPI ROOT    DIR   : ${OPENMPI_ROOT_DIR}")
MESSAGE("              OPENMPI INCLUDE DIR   : ${OPENMPI_INCLUDE_DIR}")
MESSAGE("              OPENMPI LIBRARY DIR   : ${OPENMPI_LIBRARY_DIR}")
MESSAGE("              OPENMPI LIBS          : ${HANDLE_OPENMPI_LIBS}")
MESSAGE("              OPENMPI COMPILER FLAGS: ${OPENMPI_COMPILER_FLAGS}")


# ---- Assert package handles -----------------------------------------------------------------------#
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OPENMPI "> MISSING OPENMPI LIB HANDLES." HANDLE_OPENMPI_LIBS)

# ---- Link libraries -------------------------------------------------------------------------------#
LINK_LIBRARIES(${HANDLE_OPENMPI_LIBS})