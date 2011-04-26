file(GLOB cpp_files *.cpp)
file(GLOB hpp_files *.hpp)
add_library(scanmatch SHARED ${cpp_files})

target_link_libraries(scanmatch ${IPP_LIBS})

# set the library API version.  Increment this every time the public API
# changes.
set_target_properties(scanmatch PROPERTIES SOVERSION 1)

pods_use_pkg_config_packages(scanmatch bot2-core)

# make the shared library public
pods_install_libraries(scanmatch)
#make the headers public
pods_install_headers(${hpp_files} DESTINATION scanmatch)

pods_install_pkg_config_file(scanmatch
    LIBS -lscanmatch -lpthread -l${LCMTYPES_LIBS} -lz
    CFLAGS 
    REQUIRES bot2-core bot2-lcmgl-client
    VERSION 0.0.1)
