SET(SUBDIRPATH material/fluid/linear) 
SET(OPTION_LABEL BUILD_Material_Fluid/linear)

SET(CURRENT_DIR ${SOURCE_ELPASO}/${SUBDIRPATH})

OPTION(${OPTION_LABEL} "${CURRENT_DIR}" ON)
IF(${OPTION_LABEL})
   COLLECT_PACKAGE_DATA( ${CURRENT_DIR} ${SUBDIRPATH} ALL_SOURCES)
ENDIF(${OPTION_LABEL})