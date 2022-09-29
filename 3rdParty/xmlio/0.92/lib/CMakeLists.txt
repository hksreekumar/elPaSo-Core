#################################################################
###   XMLIO                                                   ###
#################################################################
###   Marco Schauer
###   22.10.2008

MESSAGE("> Build XMLIO ...")

cmake_minimum_required(VERSION 2.4)

#################################################################
###   ROOT-Directory                                          ###
#################################################################

  SET(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../..)

#################################################################
###   PROJECT                                                 ###
#################################################################

#################################################################
###   WINDOWS                                                 ###
#################################################################
IF(INFAM_FLAG_WINDOWS)

  IF(INFAM_FLAG_32BIT)
    PROJECT (xmlio)
    SET(XMLIO_CXX_FLAGS_DEBUG   "${XMLIO_CXX_FLAGS_DEBUG}   /MT /Gm /Wp64 /ZI" CACHE STRING   "XMLIO_CXX_FLAGS_DEBUG"   FORCE)
    SET(XMLIO_CXX_FLAGS_RELEASE "${XMLIO_CXX_FLAGS_RELEASE} /MT /Wp64 /O2 /LTCG" CACHE STRING "XMLIO_CXX_FLAGS_RELEASE" FORCE)

    SET(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   ${XMLIO_CXX_FLAGS_DEBUG}")
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${XMLIO_CXX_FLAGS_RELEASE}")
  ENDIF(INFAM_FLAG_32BIT)

  IF(INFAM_FLAG_64BIT)
    PROJECT (xmlio_x64)
    SET(XMLIO_CXX_FLAGS_DEBUG   "${XMLIO_CXX_FLAGS_DEBUG}   /MT /Gm /Wp64 /ZI" CACHE STRING   "XMLIO_CXX_FLAGS_DEBUG"   FORCE)
    SET(XMLIO_CXX_FLAGS_RELEASE "${XMLIO_CXX_FLAGS_RELEASE} /MT /Wp64 /O2 /LTCG" CACHE STRING "XMLIO_CXX_FLAGS_RELEASE" FORCE)

    SET(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   ${XMLIO_CXX_FLAGS_DEBUG}")
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${XMLIO_CXX_FLAGS_RELEASE}")
  ENDIF(INFAM_FLAG_64BIT)

  SET(LIBRARY_NAME xmlio)     
#################################################################
###   MACORS                                                  ###
#################################################################


#################################################################
###   OUTPUT_PATH                                             ###
#################################################################

  IF(INFAM_FLAG_32BIT)
    SET(LIBRARY_OUTPUT_PATH ${ELPASO_SOURCE_DIR}/3rdParty/xmlio/0.92/lib)
  ENDIF(INFAM_FLAG_32BIT)
  IF(ELPASO_OS_LIN_x86_64)
    SET(LIBRARY_OUTPUT_PATH ${ELPASO_SOURCE_DIR}/3rdParty/xmlio/0.92/lib)
  ENDIF(ELPASO_OS_LIN_x86_64)

#################################################################
###   PACKAGES                                                ###
#################################################################
  
  INCLUDE(${ELPASO_SOURCE_DIR}/3rdParty/xmlio/0.92/include/CMakePackage.cmake)
  SET(XMLIO_INCLUDE_DIR ${ELPASO_SOURCE_DIR}/3rdParty/xmlio/0.92/include)
  INCLUDE_DIRECTORIES(${XMLIO_INCLUDE_DIR})

  INCLUDE(${ELPASO_SOURCE_DIR}/3rdParty/xmlio/0.92/src/CMakePackage.cmake)

#################################################################
###   MAIN                                                    ###
#################################################################

  #SET(SPECIFIC_FILES ${SOURCE_3RDPARTY}/xmlio/0.92/src/gzstream.cpp)
  #SET(ALL_SOURCES ${ALL_SOURCES}
  #                ${SPECIFIC_FILES}
  #   )
  #SOURCE_GROUP(main FILES ${SPECIFIC_FILES})

#################################################################
###   DEFINITIONS                                             ###
#################################################################

  IF(INFAM_FLAG_32BIT)
    ADD_DEFINITIONS(-D_LIB -D_CONSOLE -D_MT)
  ENDIF(INFAM_FLAG_32BIT)
  IF(ELPASO_OS_LIN_x86_64)
    ADD_DEFINITIONS(-D_LIB -D_CONSOLE -D_MT)
  ENDIF(ELPASO_OS_LIN_x86_64)

#################################################################
###   ADDITIONAL_MAKE_CLEAN_FILES                             ###
#################################################################

  SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${GENERATED_FILES}")
  
#################################################################
###   EXCECUTABLE                                             ###
#################################################################

  ADD_LIBRARY( ${LIBRARY_NAME} 
               ${ALL_SOURCES}
             )

ELSE(INFAM_FLAG_WINDOWS)
  #MESSAGE(FATAL_ERROR "windows not found!")
ENDIF(INFAM_FLAG_WINDOWS)

#################################################################
###   LINUX                                                   ###
#################################################################
IF(ELPASO_OS_LIN_x86_64)
  SET(XMLIO_CXX_FLAGS_DEBUG   "${XMLIO_CXX_FLAGS_DEBUG}   /MT /Gm /Wp64 /ZI" CACHE STRING   "XMLIO_CXX_FLAGS_DEBUG"   FORCE)
  SET(XMLIO_CXX_FLAGS_RELEASE "${XMLIO_CXX_FLAGS_RELEASE} /MT /Wp64 /O2 /LTCG" CACHE STRING "XMLIO_CXX_FLAGS_RELEASE" FORCE)

  SET(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   ${XMLIO_CXX_FLAGS_DEBUG}")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${XMLIO_CXX_FLAGS_RELEASE}")

  SET(LIBRARY_NAME xmlio)     
#################################################################
###   MACORS                                                  ###
#################################################################
INCLUDE(${ELPASO_SOURCE_DIR}/cmake/macro/CMakeCollectPackageData.cmake)

#################################################################
###   OUTPUT_PATH                                             ###
#################################################################

  SET(LIBRARY_OUTPUT_PATH ${ELPASO_SOURCE_DIR}/3rdParty/xmlio/0.92/lib)

#################################################################
###   PACKAGES                                                ###
#################################################################
  
  INCLUDE(${ELPASO_SOURCE_DIR}/3rdParty/xmlio/0.92/include/CMakePackage.cmake)
  SET(XMLIO_INCLUDE_DIR ${ELPASO_SOURCE_DIR}/3rdParty/xmlio/0.92/include)
  INCLUDE_DIRECTORIES(${XMLIO_INCLUDE_DIR})

  INCLUDE(${ELPASO_SOURCE_DIR}/3rdParty/xmlio/0.92/src/CMakePackage.cmake)

#################################################################
###   DEFINITIONS                                             ###
#################################################################

  IF(INFAM_FLAG_32BIT)
    ADD_DEFINITIONS(-D_LIB -D_CONSOLE -D_MT)
  ENDIF(INFAM_FLAG_32BIT)
  IF(ELPASO_OS_LIN_x86_64)
    ADD_DEFINITIONS(-D_LIB -D_CONSOLE -D_MT)
  ENDIF(ELPASO_OS_LIN_x86_64)

#################################################################
###   ADDITIONAL_MAKE_CLEAN_FILES                             ###
#################################################################

  SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${GENERATED_FILES}")
  
#################################################################
###   EXCECUTABLE                                             ###
#################################################################

  ADD_LIBRARY( ${LIBRARY_NAME} 
               ${ALL_SOURCES}
             )

ELSE(ELPASO_OS_LIN_x86_64)
  #MESSAGE(FATAL_ERROR "linux not found!")
ENDIF(ELPASO_OS_LIN_x86_64)
