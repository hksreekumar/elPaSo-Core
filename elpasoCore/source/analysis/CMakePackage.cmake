#include analysis
SET(SUBDIRPATH analysis) 
SET(OPTION_LABEL BUILD_Analysis)

SET(CURRENT_DIR ${SOURCE_ELPASO}/${SUBDIRPATH})

OPTION(${OPTION_LABEL} "${CURRENT_DIR}" ON)
IF(${OPTION_LABEL})
   COLLECT_PACKAGE_DATA( ${CURRENT_DIR} ${SUBDIRPATH} ALL_SOURCES)
ENDIF(${OPTION_LABEL})

#additional includes
  #INCLUDE(${SOURCE_ELPASO}/analysis/eigenvalue/CMakePackage.cmake)
  INCLUDE(${SOURCE_ELPASO}/analysis/frequency/CMakePackage.cmake)
  #INCLUDE(${SOURCE_ELPASO}/analysis/geoopt/CMakePackage.cmake)
  #INCLUDE(${SOURCE_ELPASO}/analysis/static/CMakePackage.cmake)
  #INCLUDE(${SOURCE_ELPASO}/analysis/time/CMakePackage.cmake)
  #INCLUDE(${SOURCE_ELPASO}/analysis/mor/CMakePackage.cmake)