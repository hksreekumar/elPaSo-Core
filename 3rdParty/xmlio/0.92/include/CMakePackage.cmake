SET(SUBDIRPATH xmlio/0.92/include)
SET(OPTION_LABEL BUILD_xmlio_include)

SET(CURRENT_DIR ${ELPASO_SOURCE_DIR}/3rdParty/${SUBDIRPATH})

OPTION(${OPTION_LABEL} "${CURRENT_DIR}" ON)
IF(${OPTION_LABEL})
   COLLECT_PACKAGE_DATA( ${CURRENT_DIR} ${SUBDIRPATH} ALL_SOURCES)
ENDIF(${OPTION_LABEL})