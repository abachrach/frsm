file(GLOB cpp_files *.cpp)
file(GLOB hpp_files *.h*)
add_library(frsm SHARED ${cpp_files})

set(REQUIRED_LIBS ${IPP_LIBS} -lz)

set(REQUIRED_PKGS glib-2.0)
find_package(PkgConfig REQUIRED)
pkg_check_modules(LCM lcm)
if(NOT LCM_FOUND)
    message("LCM not found.  Disabling LCM support")
    add_definitions(-DNO_LCM)    
else()
    set(REQUIRED_PKGS ${REQUIRED_PKGS} lcm lcmtypes_frsm)
endif()

pkg_check_modules(BOT_LCMGL bot2-lcmgl-client)
if(NOT BOT_LCMGL_FOUND)
    message("bot-lcmgl not found.  Disabling drawing support")
    add_definitions(-DNO_BOT_LCMGL)    
else()
    set(REQUIRED_PKGS ${REQUIRED_PKGS} bot2-lcmgl-client)
endif()


target_link_libraries(frsm ${REQUIRED_LIBS})
pods_use_pkg_config_packages(frsm ${REQUIRED_PKGS})

# make the shared library public
pods_install_libraries(frsm)
#make the headers public
pods_install_headers(${hpp_files} DESTINATION frsm)

pods_install_pkg_config_file(frsm
    LIBS -lfrsm ${REQUIRED_LIBS}
    CFLAGS 
    REQUIRES ${REQUIRED_PKGS} 
    VERSION 0.0.1)

add_subdirectory(test)
