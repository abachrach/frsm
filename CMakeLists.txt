cmake_minimum_required(VERSION 2.6.0)

set(POD_NAME frsm)
include(cmake/pods.cmake)

find_package(PkgConfig REQUIRED)
pkg_check_modules(LCM lcm)
if(NOT LCM_FOUND)
    message("LCM not found.  Disabling LCM support")
    add_definitions(-DNO_LCM)
else()
    # automatically build LCM types.  This also defines a number of CMake
    # variables, see cmake/lcmtypes.cmake for details
    include(cmake/lcmtypes.cmake)
    lcmtypes_build()
endif()

#include(cmake/ipp.cmake)

add_subdirectory(src/libfrsm)
if(LCM_FOUND)
    add_subdirectory(src/frsm-lcm)
endif()