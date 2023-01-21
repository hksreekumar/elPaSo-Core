SET(SUBDIRPATH misc) 
SET(OPTION_LABEL BUILD_Miscellaneous)

SET(CURRENT_DIR ${SOURCE_ELPASO}/${SUBDIRPATH})

OPTION(${OPTION_LABEL} "${CURRENT_DIR}" ON)
IF(${OPTION_LABEL})
   COLLECT_PACKAGE_DATA( ${CURRENT_DIR} ${SUBDIRPATH} ALL_SOURCES)
ENDIF(${OPTION_LABEL})

#additional includes
  #INCLUDE(${SOURCE_ELPASO}/misc/interface/CMakePackage.cmake)
  INCLUDE(${SOURCE_ELPASO}/misc/log/CMakePackage.cmake)
  INCLUDE(${SOURCE_ELPASO}/misc/nodalforce/CMakePackage.cmake)
  INCLUDE(${SOURCE_ELPASO}/misc/nodalmoment/CMakePackage.cmake)
  INCLUDE(${SOURCE_ELPASO}/misc/nodalpressure/CMakePackage.cmake)
  INCLUDE(${SOURCE_ELPASO}/misc/postprocess/CMakePackage.cmake)
  INCLUDE(${SOURCE_ELPASO}/misc/nodalvalues/CMakePackage.cmake)
  INCLUDE(${SOURCE_ELPASO}/misc/parser/CMakePackage.cmake)
  INCLUDE(${SOURCE_ELPASO}/misc/hdf5/CMakePackage.cmake)