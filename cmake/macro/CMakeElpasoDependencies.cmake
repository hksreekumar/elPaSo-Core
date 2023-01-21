#---------------------------------------------------------------------------#
# elPaSo - CMake Project for Module elpasoC
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut f�r Akustik, Technische Universit�t Braunschweig
#
# CMake project rewritten from THE INFAM PROJECT, Marco Schauer, 08.03.2011
#---------------------------------------------------------------------------#

MACRO(SET_ELPASO_DEPENDENCIES)
  SITE_NAME(MACHINE_NAME)
  MESSAGE("> Set SET_ELPASO_DEPENDENCIES on ${MACHINE_NAME} ...")
  INCLUDE(./cmake/${MACHINE_NAME}.config.cmake)
  DEPENDENCY_CONFIG()
ENDMACRO(SET_ELPASO_DEPENDENCIES)
