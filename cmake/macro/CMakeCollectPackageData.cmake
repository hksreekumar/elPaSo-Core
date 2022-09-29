#---------------------------------------------------------------------------#
# elPaSo - CMake Project
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut für Akustik, Technische Universität Braunschweig
#
# CMake project rewritten from THE INFAM PROJECT, Marco Schauer, 22.10.2008
#---------------------------------------------------------------------------#

#################################################################################
###                                                                           ###
###                              WITH SUBFOLDERS                              ###
###                                                                           ###
#################################################################################

IF(NOT WITH_SUBFOLDERS)
   SET(WITH_SUBFOLDERS NO)
ENDIF(NOT WITH_SUBFOLDERS)


#################################################################################
### COLLECT_PACKAGE_DATA( currentDir  sourceGroupWithSubfolders outFiles)     ###
### collects header and cpp file of current dir and add them to "outfiles"    ###
### all files will be put to the SOURCE_GROUP-folder "sourceGroupName"        ###
### and this one will be with subfolders if  WITH_SUBFOLDERS==YES             ###
#################################################################################
MACRO(COLLECT_PACKAGE_DATA currentDir sourceGroupName outFiles)
  FILE(GLOB HEADER_FILES ${currentDir}/*.h   )
  FILE(GLOB C_FILES      ${currentDir}/*.c   )
  FILE(GLOB CPP_FILES    ${currentDir}/*.cpp )
  FILE(GLOB HPP_FILES    ${currentDir}/*.hpp )
  FILE(GLOB F_FILES      ${currentDir}/*.f   )
  FILE(GLOB FD_FILES     ${currentDir}/*.fd  )
  FILE(GLOB F90_FILES    ${currentDir}/*.f90 )
  FILE(GLOB PCL_FILES    ${currentDir}/*.pcl )

  IF(INFAM_PACKAGE_DEFINTIONS)
    SET_SOURCE_FILES_PROPERTIES(${CPP_FILES} PROPERTIES COMPILE_FLAGS ${INFAM_PACKAGE_DEFINTIONS})
  ENDIF(INFAM_PACKAGE_DEFINTIONS)
  
  SET(tmpSourceGroupName ${sourceGroupName})
  IF(${WITH_SUBFOLDERS} MATCHES YES)  
    STRING(REGEX REPLACE "/" "\\\\" tmpSourceGroupName ${sourceGroupName})
  ENDIF(${WITH_SUBFOLDERS} MATCHES YES)
  SOURCE_GROUP(${tmpSourceGroupName} FILES ${HEADER_FILES} ${C_FILES} ${CPP_FILES} ${HPP_FILES} ${F_FILES} ${FD_FILES} ${F90_FILES} ${PCL_FILES})
  #SET(${outFiles} ${${outFiles}} ${HEADER_FILES} ${C_FILES} ${CPP_FILES} ${HPP_FILES} ${F_FILES} ${FD_FILES} ${F90_FILES} ${PCL_FILES})
  LIST(APPEND ${outFiles} ${HEADER_FILES} ${C_FILES} ${CPP_FILES} ${HPP_FILES} ${F_FILES} ${FD_FILES} ${F90_FILES} ${PCL_FILES} )
ENDMACRO(COLLECT_PACKAGE_DATA  currentDir sourceGroupName sourceGroupWithSubfolders outFiles)
