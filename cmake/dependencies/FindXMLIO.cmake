#---------------------------------------------------------------------------#
# elPaSo - CMake Project for XMLIO
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut für Akustik, Technische Universität Braunschweig
#---------------------------------------------------------------------------#

# ---- Start point for XMLIO libs finding -----------------------------------------------------------#
MESSAGE("> Finding XMLIO ...")

IF(${ELPASO_OS_LIN_x86_64} MATCHES FOUND)
    # ---- Set Root Directory -----------------------------------------------------------------------#
    SET(XMLIO_ROOT_DIR "${ELPASO_SOURCE_DIR}/3rdParty/xmlio/${XMLIO_VERSION}")             # XMLIO_DIR
    
    IF(EXISTS ${XMLIO_ROOT_DIR})
        # ---- Set each library ---------------------------------------------------------------------#
        SET(XMLIO_INCLUDE_DIR "${XMLIO_ROOT_DIR}/include")    # XMLIO_INCLUDE_DIR
        SET(XMLIO_LIBRARY_DIR "${XMLIO_ROOT_DIR}/lib")        # XMLIO_LIBRARY_DIR
        
        # ---- Link and include directories ---------------------------------------------------------#
        LINK_DIRECTORIES(${XMLIO_LIBRARY_DIR})
        INCLUDE_DIRECTORIES(${XMLIO_INCLUDE_DIR})

        # ---- Find each Library --------------------------------------------------------------------#
        SET(HANDLE_XMLIO_LIBS libxmlio.a)

        # ---- Set FLAGS for XMLIO ------------------------------------------------------------------#
        SET(XMLIO_COMPILER_FLAGS "")
    ELSE()
        MESSAGE(FATAL_ERROR "> XMLIO_ROOT_DIR not found.")
    ENDIF()

ENDIF()

# ---- Display summary ------------------------------------------------------------------------------#
MESSAGE("              XMLIO ROOT    DIR   : ${XMLIO_ROOT_DIR}")
MESSAGE("              XMLIO INCLUDE DIR   : ${XMLIO_INCLUDE_DIR}")
MESSAGE("              XMLIO LIBRARY DIR   : ${XMLIO_LIBRARY_DIR}")
MESSAGE("              XMLIO LIBS          : ${HANDLE_XMLIO_LIBS}")
MESSAGE("              XMLIO COMPILER FLAGS: ${XMLIO_COMPILER_FLAGS}")

# ---- Link libraries -------------------------------------------------------------------------------#
LINK_LIBRARIES(${HANDLE_XMLIO_LIBS})