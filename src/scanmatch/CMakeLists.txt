file(GLOB cpp_files *.cpp)
file(GLOB hpp_files *.hpp)
add_library(scanmatch SHARED ${cpp_files})

target_link_libraries(scanmatch ${IPP_LIBS})

# set the library API version.  Increment this every time the public API
# changes.
set_target_properties(scanmatch PROPERTIES SOVERSION 1)

# make the shared library public
pods_install_libraries(scanmatch)
#make the headers public
pods_install_headers(${hpp_files} DESTINATION scanmatch)

pods_install_pkg_config_file(scanmatch
    LIBS -lscanmatch -lpthread -l${LCMTYPES_LIBS}
    CFLAGS 
    REQUIRES opencv
    VERSION 0.0.1)
