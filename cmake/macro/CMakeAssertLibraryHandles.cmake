#---------------------------------------------------------------------------#
# elPaSo - CMake Project
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut für Akustik, Technische Universität Braunschweig
#
# CMake project rewritten from THE INFAM PROJECT, Marco Schauer, 08.03.2011
#---------------------------------------------------------------------------#

MACRO(ASSERT_LIB_HANDLES ${HANDLE_LIBS})
    FOREACH(HANDLE_LIB ${HANDLE_LIBS})
        # ---- Check if the library exists ------------------------------------------------------#
        IF(NOT HANDLE_LIB)
            MESSAGE(FATAL_ERROR "> LIB NOT FOUND!")
        ENDIF()
    ENDFOREACH(HANDLE_LIB)
ENDMACRO(ASSERT_LIB_HANDLES ${HANDLE_LIBS})