#---------------------------------------------------------------------------#
# elPaSo - CMake Project for GTEST
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut f�r Akustik, Technische Universit�t Braunschweig
#---------------------------------------------------------------------------#

# ---- Start point for GTEST libs finding -----------------------------------------------------------#
MESSAGE("> Finding GTEST ...")

IF(${ELPASO_OS_LIN_x86_64} MATCHES FOUND)
    # ---- Set Root Directory -----------------------------------------------------------------------#
    SET(GTEST_ROOT_DIR "${ELPASO_SOURCE_DIR}/3rdParty/googletest")          # GTEST_DIR
    
    IF(EXISTS ${GTEST_ROOT_DIR})
        # ---- Set each library ---------------------------------------------------------------------#
        SET(GTEST_INCLUDE_DIR "${GTEST_ROOT_DIR}/googletest/include;${GTEST_ROOT_DIR}/googlemock/include")    # GTEST_INCLUDE_DIR
        SET(GTEST_LIBRARY_DIR "${GTEST_ROOT_DIR}/../../project/lib")    # GTEST_INCLUDE_DIR
        
        # ---- Link and include directories ---------------------------------------------------------#
        INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIR})
        LINK_DIRECTORIES(${GTEST_LIBRARY_DIR})

        # ---- Find each Library --------------------------------------------------------------------#
        SET(HANDLE_GTEST_LIB_GTEST         libgtest.a)
        SET(HANDLE_GTEST_LIB_GMOCK         libgmock.a)
        
        SET(HANDLE_GTEST_LIBS ${HANDLE_INTELMPI_LIBS} ${HANDLE_GTEST_LIB_GTEST} ${HANDLE_GTEST_LIB_GMOCK})
        # ---- Set FLAGS for GTEST ------------------------------------------------------------------#
        SET(GTEST_COMPILER_FLAGS "")
    ELSE()
        MESSAGE(FATAL_ERROR "> GTEST_ROOT_DIR not found.")
    ENDIF()

ENDIF()

# ---- Display summary ------------------------------------------------------------------------------#
MESSAGE("              GTEST ROOT    DIR   : ${GTEST_ROOT_DIR}")
MESSAGE("              GTEST INCLUDE DIR   : ${GTEST_INCLUDE_DIR}")
MESSAGE("              GTEST LIBRARY DIR   : ${GTEST_LIBRARY_DIR}")
MESSAGE("              GTEST LIBS          : ${HANDLE_GTEST_LIBS}")
MESSAGE("              GTEST COMPILER FLAGS: ${GTEST_COMPILER_FLAGS}")

# ---- Link libraries -------------------------------------------------------------------------------#
LINK_LIBRARIES("${HANDLE_GTEST_LIBS};-lpthread")