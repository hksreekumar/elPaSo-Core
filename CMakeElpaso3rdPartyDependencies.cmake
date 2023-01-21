#---------------------------------------------------------------------------#
# elPaSo - CMake Project
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institute for Acoustics, Technische Universitaet Braunschweig
#
# CMake project rewritten from THE INFAM PROJECT, Marco Schauer, 22.10.2008
#---------------------------------------------------------------------------#

MACRO(SET_ELPASO_3RDPARTY_DEPENDENCIES)
	IF(${ELPASO_OS_LIN_x86_64} MATCHES FOUND)
		IF(LINK_XMLIO)
			SET(XMLIO_ROOT_DIR "${ELPASO_SOURCE_DIR}/3rdParty/xmlio/${XMLIO_VERSION}")             # XMLIO_DIR
			ADD_SUBDIRECTORY(${XMLIO_ROOT_DIR}/lib)
		ENDIF()
		IF(LINK_GZSTREAM)
			SET(GZSTREAM_ROOT_DIR "${ELPASO_SOURCE_DIR}/3rdParty/gzstream")						   # GSTREAM_DIR
			ADD_SUBDIRECTORY(${GZSTREAM_ROOT_DIR}/lib)
		ENDIF()
		IF(LINK_GTEST)
			MESSAGE("> Build GTEST ...")
			SET(GTEST_ROOT_DIR "${ELPASO_SOURCE_DIR}/3rdParty/googletest")                         # GTEST_DIR

			ADD_SUBDIRECTORY(${GTEST_ROOT_DIR} "${GTEST_ROOT_DIR}/lib")
			enable_testing()
		ENDIF()
	ENDIF()
ENDMACRO(SET_ELPASO_3RDPARTY_DEPENDENCIES)