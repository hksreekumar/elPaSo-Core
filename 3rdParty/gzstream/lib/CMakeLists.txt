#################################################################
###   GZSTREAM                                                ###
#################################################################
###   Marco Schauer
###   22.10.2008
MESSAGE("> Build GZSTREAM ...")
  cmake_minimum_required(VERSION 2.4)

#################################################################
###   ROOT-Directory                                          ###
#################################################################
SET(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
SET(SOURCE_3RDPARTY ${ELPASO_SOURCE_DIR}/3rdParty)

#################################################################
###   INFAM_FLAGS                                             ###
#################################################################

  #INCLUDE(${ROOT_DIR}/CMakeINFAMFlags.cmake)
  #SET_INFAM_FLAGS()

#################################################################
###   PROJECT                                                 ###
#################################################################

#################################################################
###   WINDOWS                                                 ###
#################################################################
IF(INFAM_FLAG_WINDOWS)

  IF(INFAM_FLAG_32BIT)
    PROJECT (gzstream)
    SET(GZSTREAM_CXX_FLAGS_DEBUG   "${GZSTREAM_CXX_FLAGS_DEBUG}   /MTd" CACHE STRING "MTd" FORCE)
    SET(GZSTREAM_CXX_FLAGS_RELEASE "${GZSTREAM_CXX_FLAGS_RELEASE} /MD"  CACHE STRING "MD"  FORCE)

    SET(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   ${GZSTREAM_CXX_FLAGS_DEBUG}")
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${GZSTREAM_CXX_FLAGS_RELEASE}")
  ENDIF(INFAM_FLAG_32BIT)

  IF(INFAM_FLAG_64BIT)
    PROJECT (gzstream_x64)
    SET(GZSTREAM_CXX_FLAGS_DEBUG   "${GZSTREAM_CXX_FLAGS_DEBUG}   /MTd" CACHE STRING "MTd" FORCE)
    SET(GZSTREAM_CXX_FLAGS_RELEASE "${GZSTREAM_CXX_FLAGS_RELEASE} /MD"  CACHE STRING "MD"  FORCE)

    SET(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   ${GZSTREAM_CXX_FLAGS_DEBUG}")
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${GZSTREAM_CXX_FLAGS_RELEASE}")
  ENDIF(INFAM_FLAG_64BIT)

  SET(LIBRARY_NAME gzstream)  
   
#################################################################
###   OUTPUT_PATH                                             ###
#################################################################

  IF(INFAM_FLAG_32BIT)
    SET(LIBRARY_OUTPUT_PATH ${SOURCE_3RDPARTY}/gzstream/lib)
  ENDIF(INFAM_FLAG_32BIT)
  IF(INFAM_FLAG_64BIT)
    SET(LIBRARY_OUTPUT_PATH ${SOURCE_3RDPARTY}/gzstream/lib64)
  ENDIF(INFAM_FLAG_64BIT)

#################################################################
###   MACORS                                                  ###
#################################################################
  
  INCLUDE(${ROOT_DIR}/CMakeINFAMMacros.cmake)

#################################################################
###   PACKAGES                                                ###
#################################################################
  
  INCLUDE(${SOURCE_3RDPARTY}/gzstream/include/CMakePackage.cmake)
  SET(GSTREAM_INCLUDE_DIR ${SOURCE_3RDPARTY}/gzstream/include)
  INCLUDE_DIRECTORIES(${GSTREAM_INCLUDE_DIR})

#################################################################
###   MAIN                                                    ###
#################################################################

  SET(SPECIFIC_FILES ${SOURCE_3RDPARTY}/gzstream/src/gzstream.cpp)
  SET(ALL_SOURCES ${ALL_SOURCES}
                  ${SPECIFIC_FILES}
     )
  SOURCE_GROUP(main FILES ${SPECIFIC_FILES})

#################################################################
###   CONTRIB                                                 ###
#################################################################

  SET_ZLIB(YES)

#################################################################
###   DEFINITIONS                                             ###
#################################################################
 
  IF(INFAM_FLAG_32BIT)
    ADD_DEFINITIONS(-D_LIB -D_CONSOLE -DINFAMVISUALSTUDIO)
  ENDIF(INFAM_FLAG_32BIT)
  IF(INFAM_FLAG_64BIT)
    ADD_DEFINITIONS(-D_LIB -D_CONSOLE -DINFAMVISUALSTUDIO)
  ENDIF(INFAM_FLAG_64BIT)

#################################################################
###   ADDITIONAL_MAKE_CLEAN_FILES                             ###
#################################################################

  SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${GENERATED_FILES}")
  
#################################################################
###   EXCECUTABLE                                             ###
#################################################################

  ADD_LIBRARY( ${LIBRARY_NAME}
               STATIC
               ${ALL_SOURCES}
             )

ELSE(INFAM_FLAG_WINDOWS)
  #MESSAGE(FATAL_ERROR "windows not found!")
ENDIF(INFAM_FLAG_WINDOWS)

#################################################################
###   LINUX                                                   ###
#################################################################
IF(ELPASO_OS_LIN_x86_64)

  PROJECT (gzstream)
  SET(GZSTREAM_CXX_FLAGS_DEBUG   "${GZSTREAM_CXX_FLAGS_DEBUG}   /MTd" CACHE STRING "MTd" FORCE)
  SET(GZSTREAM_CXX_FLAGS_RELEASE "${GZSTREAM_CXX_FLAGS_RELEASE} /MD"  CACHE STRING "MD"  FORCE)

  SET(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   ${GZSTREAM_CXX_FLAGS_DEBUG}")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${GZSTREAM_CXX_FLAGS_RELEASE}")

  SET(LIBRARY_NAME gzstream)  
   
#################################################################
###   OUTPUT_PATH                                             ###
#################################################################
  
  SET(LIBRARY_OUTPUT_PATH ${SOURCE_3RDPARTY}/gzstream/lib)
  
#################################################################
###   MACORS                                                  ###
#################################################################
  
  INCLUDE(${ELPASO_SOURCE_DIR}/cmake/macro/CMakeCollectPackageData.cmake)

#################################################################
###   PACKAGES                                                ###
#################################################################
  
  INCLUDE(${SOURCE_3RDPARTY}/gzstream/include/CMakePackage.cmake)
  SET(GSTREAM_INCLUDE_DIR ${SOURCE_3RDPARTY}/gzstream/include)
  INCLUDE_DIRECTORIES(${GSTREAM_INCLUDE_DIR})

#################################################################
###   MAIN                                                    ###
#################################################################

  SET(SPECIFIC_FILES ${SOURCE_3RDPARTY}/gzstream/src/gzstream.cpp)
  SET(ALL_SOURCES ${ALL_SOURCES}
                  ${SPECIFIC_FILES}
     )
  SOURCE_GROUP(main FILES ${SPECIFIC_FILES})

#################################################################
###   CONTRIB                                                 ###
#################################################################

 # SET_ZLIB(YES)

#################################################################
###   DEFINITIONS                                             ###
#################################################################

  ADD_DEFINITIONS(-D_LIB -D_CONSOLE)
 
#################################################################
###   ADDITIONAL_MAKE_CLEAN_FILES                             ###
#################################################################

  SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${GENERATED_FILES}")
  
#################################################################
###   EXCECUTABLE                                             ###
#################################################################

  ADD_LIBRARY( ${LIBRARY_NAME}
               STATIC
               ${ALL_SOURCES}
             )

ELSE(ELPASO_OS_LIN_x86_64)
  #MESSAGE(FATAL_ERROR "linux not found!")
ENDIF(ELPASO_OS_LIN_x86_64)
