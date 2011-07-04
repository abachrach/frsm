file(GLOB cpp_files *.cpp)
file(GLOB hpp_files *.hpp)
add_library(frsm SHARED ${cpp_files})

target_link_libraries(frsm ${IPP_LIBS})

# set the library API version.  Increment this every time the public API
# changes.
set_target_properties(frsm PROPERTIES SOVERSION 1)

# make the shared library public
pods_install_libraries(frsm)
#make the headers public
pods_install_headers(${hpp_files} DESTINATION frsm)

pods_install_pkg_config_file(frsm
    LIBS -lfrsm ${IPP_LIBS} -lpthread
    CFLAGS 
    REQUIRES lcmtypes_frsm
    VERSION 0.0.1)
