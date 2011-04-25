cmake_minimum_required(VERSION 2.6.0)

include(cmake/pods.cmake)

# automatically build LCM types.  This also defines a number of CMake
# variables, see cmake/lcmtypes.cmake for details
include(cmake/lcmtypes.cmake)
lcmtypes_build()

include(cmake/ipp.cmake)

add_subdirectory(src/scanmatch)
add_subdirectory(src/scan-matcher)
