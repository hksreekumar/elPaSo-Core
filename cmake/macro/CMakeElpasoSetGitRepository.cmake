#---------------------------------------------------------------------------#
# elPaSo - CMake Project
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institut für Akustik, Technische Universität Braunschweig
#
# CMake project rewritten from THE INFAM PROJECT, Marco Schauer, 08.03.2011
#---------------------------------------------------------------------------#

MACRO(SET_ELPASO_GIT)
MESSAGE("> Finding GIT details ...")
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

  execute_process(
    COMMAND git rev-list --all --count --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COUNT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    add_definitions(-DGIT_COUNT="\\"${GIT_COUNT}\\"")
    add_definitions(-DGIT_COMMIT_HASH="\\"${GIT_COMMIT_HASH}\\"")

    MESSAGE("              COUNT ${GIT_COUNT} | COMMIT HASH ${GIT_COMMIT_HASH}")
ENDMACRO(SET_ELPASO_GIT)
