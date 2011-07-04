include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${LCMTYPES_INCLUDE_DIRS})


add_executable(frsm scan_matcher.cpp)
target_link_libraries(frsm z)
pods_use_pkg_config_packages(frsm frsm lcm bot2-core bot2-lcmgl-client)
pods_install_executables(frsm)




