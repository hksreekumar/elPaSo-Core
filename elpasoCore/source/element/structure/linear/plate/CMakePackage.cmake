SET(SUBDIRPATH element/structure/linear/plate) 
SET(OPTION_LABEL BUILD_Element_Structure_Linear_Plate)

SET(CURRENT_DIR ${SOURCE_ELPASO}/${SUBDIRPATH})

OPTION(${OPTION_LABEL} "${CURRENT_DIR}" ON)
IF(${OPTION_LABEL})
   COLLECT_PACKAGE_DATA( ${CURRENT_DIR} ${SUBDIRPATH} ALL_SOURCES)
ENDIF(${OPTION_LABEL})