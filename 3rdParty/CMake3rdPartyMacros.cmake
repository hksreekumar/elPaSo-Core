#################################################################################
###                                                                           ###
###                               SET_484                                     ###
###                                                                           ###
#################################################################################
MACRO(SET_484 useStaticLibs)

  SET(OPTION_LABEL SET_484)

  IF(${useStaticLibs} MATCHES YES)
    SET(484_DIR ${SOURCE_3RDPARTY}/484)                      # 484_DIR

    IF(EXISTS ${484_DIR})
###32bit
      IF(INFAM_FLAG_WINDOWS)
        OPTION(${OPTION_LABEL} "${ARPACK_DIR}" ON)
        SET(484_LIB_DIR ${484_DIR}/lib)                             # 484_LIB_DIR
        LINK_DIRECTORIES(${484_LIB_DIR})
        LINK_LIBRARIES(_484_.lib)#;484.lib)
        #LINK_LIBRARIES(484.lib)
      ENDIF(INFAM_FLAG_WINDOWS)
##64bit
      IF(INFAM_FLAG_LINUX)
        OPTION(${OPTION_LABEL} "${484_DIR}" ON)
        SET(484_LIB_DIR ${484_DIR}/lib)                           # 484_LIB64_DIR
        LINK_DIRECTORIES(${484_LIB_DIR})
        #LINK_LIBRARIES(lib_484_.a)#;484.lib)
      ENDIF(INFAM_FLAG_LINUX)

    ELSE(EXISTS ${484_DIR})
      OPTION(${OPTION_LABEL} "${484_DIR}" OFF)
      MESSAGE("SET_484 - failed - no 484-dir!")
    ENDIF(EXISTS ${484_DIR})

  ELSE(${useStaticLibs} MATCHES YES)
    OPTION(${OPTION_LABEL} "${484_DIR}" OFF)
  ENDIF(${useStaticLibs} MATCHES YES)
ENDMACRO(SET_484 useStaticLibs)

#################################################################################
###                                                                           ###
###                              SET_ARPACK                                   ###
###                                                                           ###
#################################################################################
MACRO(SET_ARPACK useStaticLibs)

  SET(OPTION_LABEL SET_ARPACK)

  IF(${useStaticLibs} MATCHES YES)
    IF(INFAM_FLAG_WINDOWS)
      SET(ARPACK_DIR ${SOURCE_3RDPARTY}/ARPACK)                      # ARPACK_DIR

      IF(EXISTS ${ARPACK_DIR})
###32bit
        IF(INFAM_FLAG_32BIT)
          OPTION(${OPTION_LABEL} "${ARPACK_DIR}" ON)
          SET(ARPACK_LIB_DIR ${ARPACK_DIR}/lib)                             # ARPACK_LIB_DIR
          LINK_DIRECTORIES(${ARPACK_LIB_DIR})
          LINK_LIBRARIES(_arpack_.lib;ARPACK.lib)
        ENDIF(INFAM_FLAG_32BIT)
###64bit
        IF(INFAM_FLAG_64BIT)
          OPTION(${OPTION_LABEL} "${ARPACK_DIR}" ON)
          SET(ARPACK_LIB_DIR ${ARPACK_DIR}/lib64)                           # ARPACK_LIB64_DIR
          LINK_DIRECTORIES(${ARPACK_LIB_DIR})
          LINK_LIBRARIES(_arpack_.lib;ARPACK.lib)
        ENDIF(INFAM_FLAG_64BIT)

      ELSE(EXISTS ${ARPACK_DIR})
        OPTION(${OPTION_LABEL} "${ARPACK_DIR}" OFF)
        MESSAGE("SET_ARPACK - failed - no arpack-dir!")
      ENDIF(EXISTS ${ARPACK_DIR})
    ENDIF(INFAM_FLAG_WINDOWS)
	IF(INFAM_FLAG_LINUX)
	  LINK_DIRECTORIES(${ARPACK_DIR})
      LINK_LIBRARIES(libparpack_${PETSC_ARCH}.a
	                 libarpack_${PETSC_ARCH}.a
	                )
	ENDIF(INFAM_FLAG_LINUX)
  ELSE(${useStaticLibs} MATCHES YES)
    OPTION(${OPTION_LABEL} "${ARPACK_DIR}" OFF)
  ENDIF(${useStaticLibs} MATCHES YES)
ENDMACRO(SET_ARPACK useStaticLibs)

#################################################################################
###                                                                           ###
###                              SET_FFTW                                     ###
###                                                                           ###
#################################################################################
MACRO(SET_FFTW useStaticLibs)

  SET(OPTION_LABEL SET_FFTW)

  IF(${useStaticLibs} MATCHES YES)
    ADD_DEFINITIONS(-DHAVE_FFTW)    # Add definition
      
      IF(INFAM_FLAG_WINDOWS)
        SET(FFTW_DIR ${SOURCE_3RDPARTY}/fftw/3.1.2-win)                # FFTW_DIR
        OPTION(${OPTION_LABEL} "${FFTW_DIR}" ON)
        IF(EXISTS ${FFTW_DIR})
          SET(FFTW_LIB_DIR ${FFTW_DIR})                                       # FFTW_LIB_DIR
          LINK_DIRECTORIES(${FFTW_LIB_DIR})
          LINK_LIBRARIES(libfftw3-3.lib)
          INCLUDE_DIRECTORIES(${FFTW_DIR})				  # h-files
        ELSE(EXISTS ${FFTW_DIR})
          OPTION(${OPTION_LABEL} "${FFTW_DIR}" OFF)
          MESSAGE("SET_FFTW - failed - no fftw-dir!")
        ENDIF(EXISTS ${FFTW_DIR})
      ENDIF(INFAM_FLAG_WINDOWS)
      IF(INFAM_FLAG_LINUX)                        
        MESSAGE("FFTW version ${INFAM_FLAG_FFTW_VER}")
        MESSAGE("${FFTW_DIR}/include")
        MESSAGE("${FFTW_DIR}/lib")
                
        INCLUDE_DIRECTORIES(${FFTW_DIR}/include)
        LINK_DIRECTORIES(${FFTW_DIR}/lib/)
        LINK_LIBRARIES(${FFTW_LIBS})      
               
     ENDIF(INFAM_FLAG_LINUX)

  ELSE(${useStaticLibs} MATCHES YES)
    OPTION(${OPTION_LABEL} "${FFTW_DIR}" OFF)
  ENDIF(${useStaticLibs} MATCHES YES)
ENDMACRO(SET_FFTW useStaticLibs)

#################################################################################
###                                                                           ###
###                              SET_GSTREAM                                  ###
###                                                                           ###
#################################################################################
MACRO(SET_GSTREAM useStaticLibs)

  SET(OPTION_LABEL SET_GSTREAM)

  IF(${useStaticLibs} MATCHES YES)
    SET(GSTREAM_DIR ${SOURCE_3RDPARTY}/gzstream)                      # GSTREAM_DIR
    IF(EXISTS ${GSTREAM_DIR})
      OPTION(${OPTION_LABEL} "${GSTREAM_DIR}" ON)
      INCLUDE(${SOURCE_3RDPARTY}/gzstream/include/CMakePackage.cmake)
      SET(GSTREAM_INCLUDE_DIR ${GSTREAM_DIR}/include)                        # GSTREAM_INCLUDE_DIR
      INCLUDE_DIRECTORIES(${GSTREAM_INCLUDE_DIR})

      IF(INFAM_FLAG_WINDOWS)
###32bit
        IF(INFAM_FLAG_32BIT)
          SET(GSTREAM_LIB_DIR ${GSTREAM_DIR}/lib)                              # GSTREAM_LIB_DIR
          LINK_DIRECTORIES(${GSTREAM_LIB_DIR})
          LINK_LIBRARIES(gzstream.lib)
        ENDIF(INFAM_FLAG_32BIT)
###64bit
        IF(INFAM_FLAG_64BIT)
          SET(GSTREAM_LIB_DIR ${GSTREAM_DIR}/lib64)                            # GSTREAM_LIB64_DIR
          LINK_DIRECTORIES(${GSTREAM_LIB_DIR})
          LINK_LIBRARIES(gzstream.lib)
        ENDIF(INFAM_FLAG_64BIT)
      ENDIF(INFAM_FLAG_WINDOWS)
      IF(INFAM_FLAG_LINUX)
        SET(GSTREAM_LIB_DIR ${GSTREAM_DIR}/lib)                              # GSTREAM_LIB_DIR
        LINK_DIRECTORIES(${GSTREAM_LIB_DIR})
        LINK_LIBRARIES(libgzstream.a
                      -lz
                      )
      ENDIF(INFAM_FLAG_LINUX)

    ELSE(EXISTS ${GSTREAM_DIR})
      OPTION(${OPTION_LABEL} "${GSTREAM_DIR}" OFF)
      MESSAGE("SET_GSTREAM - failed - no gstream-dir!")
    ENDIF(EXISTS ${GSTREAM_DIR})

  ELSE(${useStaticLibs} MATCHES YES)
    OPTION(${OPTION_LABEL} "${GSTREAM_DIR}" OFF)
  ENDIF(${useStaticLibs} MATCHES YES)
ENDMACRO(SET_GSTREAM useStaticLibs)

#################################################################################
###                                                                           ###
###                              SET_PETSC                                    ###
###                                                                           ###
#################################################################################
MACRO(SET_PETSC useStaticLibs version complex)

  SET(OPTION_LABEL SET_PETSC)
  IF(${useStaticLibs} MATCHES YES)
    ADD_DEFINITIONS(-DHAVE_PETSC)
    IF(INFAM_FLAG_WINDOWS)
      SET(PETSC_DIR ${SOURCE_3RDPARTY}/PETSc/${version})             # PETSC_DIR
      IF(EXISTS ${PETSC_DIR})
        OPTION(${OPTION_LABEL} "${PETSC_DIR}" ON)
        SET(PETSC_INCLUDE_DIR ${PETSC_DIR}/include)                  # PETSC_INCLUDE_DIR
        IF(${complex} MATCHES YES)                                   # PETSC_LIB_DIR
          SET(PETSC_LIB_DIR ${PETSC_DIR}/lib-complex)
          MESSAGE("PETSc version ${version} complex")
        ELSE(${complex} MATCHES YES)
          SET(PETSC_LIB_DIR ${PETSC_DIR}/lib)
          MESSAGE("PETSc version ${version}")
        ENDIF(${complex} MATCHES YES)
        INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIR})
        INCLUDE_DIRECTORIES(${PETSC_INCLUDE_MPI_DIR})
        INCLUDE_DIRECTORIES(${PETSC_LIB_DIR}/include)
        LINK_DIRECTORIES(${PETSC_LIB_DIR})
        IF("${version}"  STREQUAL "2.3.3-p15")
          LINK_LIBRARIES(libpetsc.lib
                         libpetscdm.lib
                         libpetscksp.lib
                         libpetscmat.lib
                         libpetscsnes.lib
                         libpetscts.lib
                         libpetscvec.lib
                         libpetsccontrib.lib
                        )
        ELSEIF("${version}"  STREQUAL "3.0.0-p8")
  	  LINK_LIBRARIES(libpetsc.lib
                         libpetscdm.lib
                         libpetscksp.lib
                         libpetscmat.lib
                         libpetscsnes.lib
                         libpetscts.lib
                         libpetscvec.lib
                         libpetsccontrib.lib
                        )
        ELSEIF("${version}"  STREQUAL "3.1-p7")
          LINK_LIBRARIES(libpetsc.lib)
        ELSEIF("${version}"  STREQUAL "3.1-p8")
          LINK_LIBRARIES(libpetsc.lib)
        ELSE()
          MESSAGE("SET_PETSC - failed - unknown PETSC version!")
        ENDIF()
      ELSE(EXISTS ${PETSC_DIR})
        OPTION(${OPTION_LABEL} "${PETSC_DIR}" OFF)
        MESSAGE("SET_PETSC - failed - no petsc-dir!")
      ENDIF(EXISTS ${PETSC_DIR})
    ENDIF(INFAM_FLAG_WINDOWS)
    IF(INFAM_FLAG_LINUX)
      IF(${complex} MATCHES YES)                                   # PETSC_LIB_DIR
        SET(PETSC_ARCH ${PETSC_ARCH_COMPLEX})
      ELSE(${complex} MATCHES YES)
        SET(PETSC_ARCH ${PETSC_ARCH_REAL})
      ENDIF(${complex} MATCHES YES)
      MESSAGE("PETSc version ${version} ${PETSC_ARCH}")
      MESSAGE("${PETSC_DIR}/include")
      MESSAGE("${PETSC_DIR}/${PETSC_ARCH}/include")
      MESSAGE("${PETSC_DIR}/${PETSC_ARCH}/lib")

      INCLUDE_DIRECTORIES(${PETSC_DIR}/include)
      INCLUDE_DIRECTORIES(${PETSC_DIR}/${PETSC_ARCH}/include)
      LINK_DIRECTORIES(${PETSC_DIR}/${PETSC_ARCH}/lib)
      LINK_LIBRARIES(${PETSC_LIBS})
      
      IF(${complex} MATCHES YES)                                   # additional PETSC_LIBS
        LINK_LIBRARIES(${PETSC_LIBS_COMPLEX})
      ELSE(${complex} MATCHES YES)
        LINK_LIBRARIES(${PETSC_LIBS_REAL})
      ENDIF(${complex} MATCHES YES)

      OPTION(${OPTION_LABEL} "${PETSC_DIR}/${PETSC_ARCH}" ON)
    ENDIF(INFAM_FLAG_LINUX)
  ELSE(${useStaticLibs} MATCHES YES)
    OPTION(${OPTION_LABEL} "${PETSC_DIR}" OFF)
  ENDIF(${useStaticLibs} MATCHES YES)

ENDMACRO(SET_PETSC useStaticLibs)

#################################################################################
###                                                                           ###
###                              SET_SLEPC                              ###
###                                                                           ###
#################################################################################
MACRO(SET_SLEPC useStaticLibs version complex)

  SET(OPTION_LABEL SET_SLEPC)

  IF(${useStaticLibs} MATCHES YES)
    IF(INFAM_FLAG_WINDOWS)     
      SET(SLEPC_DIR ${SOURCE_3RDPARTY}/slepc/${version})             # SLEPC_DIR
      IF(EXISTS ${SLEPC_DIR})
        OPTION(${OPTION_LABEL} "${SLEPC_DIR}" ON)
        SET(SLEPC_INCLUDE_DIR ${SLEPC_DIR}/include)              # SLEPC_INCLUDE_DIR
        IF(${complex} MATCHES YES)                               # SLEPC_LIB_DIR
          SET(SLEPC_LIB_DIR ${SLEPC_DIR}/lib-complex)
          MESSAGE("SLEPc version ${version} complex")
        ELSE(${complex} MATCHES YES)
          SET(SLEPC_LIB_DIR ${SLEPC_DIR}/lib)
          MESSAGE("SLEPc version ${version}")
        ENDIF(${complex} MATCHES YES)
        INCLUDE_DIRECTORIES(${SLEPC_INCLUDE_DIR})
        LINK_DIRECTORIES(${SLEPC_LIB_DIR})
        LINK_LIBRARIES(libslepc.lib)
      ELSE(EXISTS ${SLEPC_DIR})
        OPTION(${OPTION_LABEL} "${SLEPC_DIR}" OFF)
        MESSAGE("SET_SLEPC - failed - no slepc-dir!")
      ENDIF(EXISTS ${SLEPC_DIR})
    ENDIF(INFAM_FLAG_WINDOWS)    
    IF(INFAM_FLAG_LINUX)
      IF(${complex} MATCHES YES)                                   # PETSC_LIB_DIR
        SET(PETSC_ARCH ${PETSC_ARCH_COMPLEX})
      ELSE(${complex} MATCHES YES)
        SET(PETSC_ARCH ${PETSC_ARCH_REAL})
      ENDIF(${complex} MATCHES YES)
      MESSAGE("SLEPc version ${version} ${PETSC_ARCH}")
      MESSAGE("${SLEPC_DIR}/include")
      MESSAGE("${SLEPC_DIR}/${PETSC_ARCH}/include")
      MESSAGE("${SLEPC_DIR}/${PETSC_ARCH}/lib")
      INCLUDE_DIRECTORIES(${SLEPC_DIR}/include)
      INCLUDE_DIRECTORIES(${SLEPC_DIR}/${PETSC_ARCH}/include)
      LINK_DIRECTORIES(${SLEPC_DIR}/${PETSC_ARCH}/lib)
      LINK_LIBRARIES(${SLEPC_LIBS})
      OPTION(${OPTION_LABEL} "${SLEPC_DIR}/${PETSC_ARCH}" ON)
	  IF(${NEED_ARPACK} MATCHES YES)
        SET_ARPACK(YES)
      ENDIF(${NEED_ARPACK} MATCHES YES)
    ENDIF(INFAM_FLAG_LINUX)
  ELSE(${useStaticLibs} MATCHES YES)
    OPTION(${OPTION_LABEL} "${SLEPC_DIR}" OFF)
  ENDIF(${useStaticLibs} MATCHES YES)

ENDMACRO(SET_SLEPC useStaticLibs)

#################################################################################
###                                                                           ###
###                              SET_XMLIO                              ###
###                                                                           ###
#################################################################################
MACRO(SET_XMLIO useStaticLibs)

  SET(OPTION_LABEL SET_XMLIO)

  IF(${useStaticLibs} MATCHES YES)
    IF(INFAM_FLAG_WINDOWS)
      SET(XMLIO_DIR ${SOURCE_3RDPARTY}/xmlio/0.92)                #XMLIO_DIR
   
      IF(EXISTS ${XMLIO_DIR})
###32bit
        IF(INFAM_FLAG_32BIT)
          OPTION(${OPTION_LABEL} "${XMLIO_DIR}" ON)
          SET(XMLIO_INCLUDE_DIR ${XMLIO_DIR}/include)                    # XMLIO_INCLUDE_DIR
          SET(XMLIO_LIB_DIR ${XMLIO_DIR}/lib)                            # XMLIO_LIB_DIR
          INCLUDE_DIRECTORIES(${XMLIO_INCLUDE_DIR})
          LINK_DIRECTORIES(${XMLIO_LIB_DIR})
          LINK_LIBRARIES(xmlio.lib)
        ENDIF(INFAM_FLAG_32BIT)
###64bit
        IF(INFAM_FLAG_64BIT)
          OPTION(${OPTION_LABEL} "${XMLIO_DIR}" ON)
          SET(XMLIO_INCLUDE_DIR ${XMLIO_DIR}/include)                    # XMLIO_INCLUDE_DIR
          SET(XMLIO_LIB_DIR ${XMLIO_DIR}/lib64)                          # XMLIO_LIB64_DIR
          INCLUDE_DIRECTORIES(${XMLIO_INCLUDE_DIR})
          LINK_DIRECTORIES(${XMLIO_LIB_DIR})
          LINK_LIBRARIES(xmlio.lib)
        ENDIF(INFAM_FLAG_64BIT)
      ELSE(EXISTS ${XMLIO_DIR})
        OPTION(${OPTION_LABEL} "${XMLIO_DIR}" OFF)
        MESSAGE("SET_XMLIO - failed - no xmlio-dir!")
      ENDIF(EXISTS ${XMLIO_DIR})
    ENDIF(INFAM_FLAG_WINDOWS)
    IF(INFAM_FLAG_LINUX)
      SET(XMLIO_DIR ${SOURCE_3RDPARTY}/xmlio/0.92)                #XMLIO_DIR
      SET(XMLIO_INCLUDE_DIR ${XMLIO_DIR}/include)                    # XMLIO_INCLUDE_DIR
      SET(XMLIO_LIB_DIR ${XMLIO_DIR}/lib)                            # XMLIO_LIB_DIR
      INCLUDE_DIRECTORIES(${XMLIO_INCLUDE_DIR})
      LINK_DIRECTORIES(${XMLIO_LIB_DIR})
      LINK_LIBRARIES(libxmlio.a)
    ENDIF(INFAM_FLAG_LINUX)
  ELSE(${useStaticLibs} MATCHES YES)
    OPTION(${OPTION_LABEL} "${XMLIO_DIR}" OFF)
  ENDIF(${useStaticLibs} MATCHES YES)

ENDMACRO(SET_XMLIO useStaticLibs)

#################################################################################
###                                                                           ###
###                              SET_ZLIB                               ###
###                                                                           ###
#################################################################################
MACRO(SET_ZLIB useStaticLibs)

  SET(OPTION_LABEL SET_ZLIB)

  IF(${useStaticLibs} MATCHES YES)
    SET(ZLIB_DIR ${SOURCE_3RDPARTY}/zlib123-dll)                # ZLIB_DIR
    OPTION(${OPTION_LABEL} "${ZLIB_DIR}" ON)

    IF(EXISTS ${ZLIB_DIR})
      SET(ZLIB_LIB_DIR ${ZLIB_DIR})                                    # ZLIB_LIB_DIR
      SET(ZLIB_INCLUDE_DIR ${ZLIB_DIR}/include)                        # ZLIB_INCLUDE_DIR
      SET(ZLIB_LIB_DIR ${ZLIB_DIR}/lib)                                # ZLIB_LIB_DIR

      INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
      LINK_DIRECTORIES(${ZLIB_LIB_DIR})
      LINK_LIBRARIES(zdll.lib)

    ELSE(EXISTS ${ZLIB_DIR})
      OPTION(${OPTION_LABEL} "${FFTW_DIR}" OFF)
      MESSAGE("SET_ZLIB - failed - no zlib-dir!")
    ENDIF(EXISTS ${ZLIB_DIR})
  ELSE(${useStaticLibs} MATCHES YES)
    OPTION(${OPTION_LABEL} "${ZLIB_DIR}" OFF)
  ENDIF(${useStaticLibs} MATCHES YES)

ENDMACRO(SET_ZLIB useStaticLibs)


#################################################################################
###                                                                           ###
###                              SET_HDF5                                     ###
###                                                                           ###
#################################################################################
MACRO(SET_HDF5 useStaticLibs version)
  SET(OPTION_LABEL SET_HDF5)

  IF(${useStaticLibs} MATCHES YES)
    IF(INFAM_FLAG_WINDOWS)     
      MESSAGE("HDF5 Support for Windows is not implemented!")
    ENDIF(INFAM_FLAG_WINDOWS)    
    IF(INFAM_FLAG_LINUX)
      ADD_DEFINITIONS(-DHAVE_HDF5)
      MESSAGE("HDF5 version ${version}")
      MESSAGE("${HDF5_DIR}/include")
      MESSAGE("${HDF5_DIR}/lib")
      INCLUDE_DIRECTORIES(${HDF5_DIR}/include)
      INCLUDE_DIRECTORIES(${HDF5_DIR}/lib)
      LINK_DIRECTORIES(${HDF5_DIR}/lib)
      LINK_LIBRARIES(${HDF5_LIBS})
      OPTION(${OPTION_LABEL} "${HDF5_DIR}" ON)
    ENDIF(INFAM_FLAG_LINUX)
  ELSE(${useStaticLibs} MATCHES YES)
    OPTION(${OPTION_LABEL} "${HDF5_DIR}" OFF)
  ENDIF(${useStaticLibs} MATCHES YES)

ENDMACRO(SET_HDF5 useStaticLibs)
