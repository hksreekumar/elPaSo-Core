#---------------------------------------------------------------------------#
# elPaSo - CMake Project for GZSTREAM
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut f�r Akustik, Technische Universit�t Braunschweig
#---------------------------------------------------------------------------#

# ---- Start point for GZSTREAM libs finding -----------------------------------------------------------#
MESSAGE("> Finding GZSTREAM ...")

IF(${ELPASO_OS_LIN_x86_64} MATCHES FOUND)
    # ---- Set Root Directory -----------------------------------------------------------------------#
    SET(GZSTREAM_ROOT_DIR "${ELPASO_SOURCE_DIR}/3rdParty/gzstream")             # GZSTREAM_DIR
    
    IF(EXISTS ${GZSTREAM_ROOT_DIR})
        # ---- Set each library ---------------------------------------------------------------------#
        SET(GZSTREAM_INCLUDE_DIR "${GZSTREAM_ROOT_DIR}/include")    # GZSTREAM_INCLUDE_DIR
        SET(GZSTREAM_LIBRARY_DIR "${GZSTREAM_ROOT_DIR}/lib")        # GZSTREAM_LIBRARY_DIR
        
        # ---- Link and include directories ---------------------------------------------------------#
        LINK_DIRECTORIES(${GZSTREAM_LIBRARY_DIR})
        INCLUDE_DIRECTORIES(${GZSTREAM_INCLUDE_DIR})

        # ---- Find each Library --------------------------------------------------------------------#
        SET(HANDLE_GZSTREAM_LIBS libgzstream.a)

        # ---- Set FLAGS for GZSTREAM ------------------------------------------------------------------#
        SET(GZSTREAM_COMPILER_FLAGS "")
    ELSE()
        MESSAGE(FATAL_ERROR "> GZSTREAM_ROOT_DIR not found.")
    ENDIF()

ENDIF()

# ---- Display summary ------------------------------------------------------------------------------#
MESSAGE("              GZSTREAM ROOT    DIR   : ${GZSTREAM_ROOT_DIR}")
MESSAGE("              GZSTREAM INCLUDE DIR   : ${GZSTREAM_INCLUDE_DIR}")
MESSAGE("              GZSTREAM LIBRARY DIR   : ${GZSTREAM_LIBRARY_DIR}")
MESSAGE("              GZSTREAM LIBS          : ${HANDLE_GZSTREAM_LIBS}")
MESSAGE("              GZSTREAM COMPILER FLAGS: ${GZSTREAM_COMPILER_FLAGS}")

# ---- Link libraries -------------------------------------------------------------------------------#
LINK_LIBRARIES(${HANDLE_GZSTREAM_LIBS} -lz)

# ---- Link flags -------------------------------------------------------------------------------#
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GZSTREAM_COMPILER_FLAGS}")