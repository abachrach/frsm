include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${LCMTYPES_INCLUDE_DIRS})

#download the concise-args header
set(CONCISE_ARGS_FILE ${CMAKE_CURRENT_SOURCE_DIR}/ConciseArgs.hpp)
if(NOT EXISTS ${CONCISE_ARGS_FILE})
    file(DOWNLOAD "http://concise-args.googlecode.com/files/ConciseArgs.hpp" ${CONCISE_ARGS_FILE} SHOW_PROGRESS)
endif()


add_executable(frsm scan_matcher.cpp)
target_link_libraries(frsm z)
pods_use_pkg_config_packages(frsm frsm lcm bot2-core bot2-lcmgl-client)
pods_install_executables(frsm)




