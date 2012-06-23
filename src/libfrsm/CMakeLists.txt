file(GLOB cpp_files *.cpp)
file(GLOB hpp_files *.hpp)
add_library(frsm SHARED ${cpp_files})

target_link_libraries(frsm ${IPP_LIBS} -lpthread)

set(REQUIRED_LIBS lcm bot2-lcmgl-client occ-map lcmtypes_frsm)
pods_use_pkg_config_packages(frsm ${REQUIRED_LIBS})

# make the shared library public
pods_install_libraries(frsm)
#make the headers public
pods_install_headers(${hpp_files} DESTINATION frsm)

pods_install_pkg_config_file(frsm
    LIBS -lfrsm ${IPP_LIBS} -lpthread
    CFLAGS 
    REQUIRES ${REQUIRED_LIBS} 
    VERSION 0.0.1)
