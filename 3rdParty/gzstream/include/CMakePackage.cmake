SET(SUBDIRPATH gzstream/include)
SET(OPTION_LABEL BUILD_gzstream)

SET(CURRENT_DIR ${SOURCE_3RDPARTY}/${SUBDIRPATH})

OPTION(${OPTION_LABEL} "${CURRENT_DIR}" ON)
IF(${OPTION_LABEL})
   COLLECT_PACKAGE_DATA( ${CURRENT_DIR} ${SUBDIRPATH} ALL_SOURCES)
ENDIF(${OPTION_LABEL})