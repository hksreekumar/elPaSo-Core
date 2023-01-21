#################################################################
###   ALL INPORTANT FRINK PACKAGES                            ###
#################################################################

SET(WITH_SUBFOLDERS YES)

# basics
INCLUDE(${SOURCE_ELPASO}/basics/CMakePackage.cmake)

# analysis
IF(INCLUDE_ANALYSIS_ESSENTIALS)
  INCLUDE(${SOURCE_ELPASO}/analysis/CMakePackage.cmake)  
ENDIF()

# boundarycondition
IF(INCLUDE_BC_ESSENTIALS)
  INCLUDE(${SOURCE_ELPASO}/bc/CMakePackage.cmake)
ENDIF()

# elements
IF(INCLUDE_ELEMENT_ESSENTIALS)
  INCLUDE(${SOURCE_ELPASO}/element/CMakePackage.cmake)
  ENDIF()

# material
IF(INCLUDE_MATERIAL_ESSENTIALS)
  INCLUDE(${SOURCE_ELPASO}/material/CMakePackage.cmake)
ENDIF()

# math
INCLUDE(${SOURCE_ELPASO}/math/CMakePackage.cmake)

# misc
INCLUDE(${SOURCE_ELPASO}/misc/CMakePackage.cmake)

# fe data
INCLUDE(${SOURCE_ELPASO}/fedata/CMakePackage.cmake)

# shape
INCLUDE(${SOURCE_ELPASO}/shape/CMakePackage.cmake)

#test
IF(ELPASO_TEST)
  INCLUDE(${SOURCE_ELPASO}/test/CMakePackage.cmake)
ENDIF()

#patran
#INCLUDE(${SOURCE_ELPASO}/../patran/CMakePackage.cmake)
