#---------------------------------------------------------------------------#
# elPaSo - CMake Project for HDF5
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut für Akustik, Technische Universität Braunschweig
#---------------------------------------------------------------------------#

# ---- Start point for HDF5 libs finding -----------------------------------------------------------#
MESSAGE("> Finding HDF5 ...")

IF(${ELPASO_OS_LIN_x86_64} MATCHES FOUND)
    # ---- Set Root Directory -----------------------------------------------------------------------#
    SET(HDF5_ROOT_DIR "${ELPASO_LIB_DIR}/hdf5-${HDF5_VERSION}/${ELPASO_COMPILER_OPT}")             # HDF5_DIR
    
    IF(EXISTS ${HDF5_ROOT_DIR})
        # ---- Set each library ---------------------------------------------------------------------#
        SET(HDF5_INCLUDE_DIR "${HDF5_ROOT_DIR}/include")    # HDF5_INCLUDE_DIR
        SET(HDF5_LIBRARY_DIR "${HDF5_ROOT_DIR}/lib")        # HDF5_LIBRARY_DIR
        
        # ---- Link and include directories ---------------------------------------------------------#
        LINK_DIRECTORIES(${HDF5_LIBRARY_DIR})
        INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})

        # ---- Find each Library --------------------------------------------------------------------#
        SET(HANDLE_HDF5_LIBS)

        #FIND_LIBRARY(HANDLE_HDF5_LIB_HDF5         libhdf5.so          ${HDF5_LIBRARY_DIR})
        #FIND_LIBRARY(HANDLE_HDF5_LIB_HDF5CPP      libhdf5_cpp.so      ${HDF5_LIBRARY_DIR})
        #FIND_LIBRARY(HANDLE_HDF5_LIB_Z            libz.so             ${HDF5_LIBRARY_DIR} NO_DEFAULT_PATH)
        #FIND_LIBRARY(HANDLE_HDF5_LIB_SZIP         libszip.so          ${HDF5_LIBRARY_DIR})
        #FIND_LIBRARY(HANDLE_HDF5_LIB_HL           libhdf5_hl.so       ${HDF5_LIBRARY_DIR})
        #FIND_LIBRARY(HANDLE_HDF5_LIB_HLCPP        libhdf5_hl_cpp.so   ${HDF5_LIBRARY_DIR})
        #FIND_LIBRARY(HANDLE_HDF5_LIB_TOOLS        libhdf5_tools.so    ${HDF5_LIBRARY_DIR})

        #SET(HANDLE_HDF5_LIBS ${HANDLE_HDF5_LIBS} ${HANDLE_HDF5_LIB_HDF5} ${HANDLE_HDF5_LIB_HDF5CPP} ${HANDLE_HDF5_LIB_Z} ${HANDLE_HDF5_LIB_SZIP} ${HANDLE_HDF5_LIB_HL} ${HANDLE_HDF5_LIB_HLCPP} ${HANDLE_HDF5_LIB_TOOLS})
        SET(HANDLE_HDF5_LIBS "")

        FIND_LIBRARY(HANDLE_HDF5_LIB_HDF5         libhdf5.so          ${HDF5_LIBRARY_DIR})
        SET(HANDLE_HDF5_LIBS ${HANDLE_HDF5_LIBS} ${HANDLE_HDF5_LIB_HDF5})

        # ---- Set FLAGS for HDF5 ------------------------------------------------------------------#
        SET(HDF5_COMPILER_FLAGS "")
    ELSE()
        MESSAGE(FATAL_ERROR "> HDF5_ROOT_DIR not found.")
    ENDIF()

ENDIF()

# ---- Display summary ------------------------------------------------------------------------------#
MESSAGE("              HDF5 VERSION       : ${HDF5_VERSION}")
MESSAGE("              HDF5 ROOT    DIR   : ${HDF5_ROOT_DIR}")
MESSAGE("              HDF5 INCLUDE DIR   : ${HDF5_INCLUDE_DIR}")
MESSAGE("              HDF5 LIBRARY DIR   : ${HDF5_LIBRARY_DIR}")
MESSAGE("              HDF5 LIBS          : ${HANDLE_HDF5_LIBS}")
MESSAGE("              HDF5 COMPILER FLAGS: ${HDF5_COMPILER_FLAGS}")

# ---- Assert package handles -----------------------------------------------------------------------#
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5 "> MISSING HDF5 LIB HANDLES." HANDLE_HDF5_LIBS)

# ---- Link libraries -------------------------------------------------------------------------------#
LINK_LIBRARIES(${HANDLE_HDF5_LIBS})