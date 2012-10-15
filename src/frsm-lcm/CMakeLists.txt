include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${LCMTYPES_INCLUDE_DIRS})

#download the concise-args header
set(CONCISE_ARGS_FILE ${CMAKE_CURRENT_SOURCE_DIR}/ConciseArgs.hpp)
if(NOT EXISTS ${CONCISE_ARGS_FILE})
    file(DOWNLOAD "http://concise-args.googlecode.com/files/ConciseArgs.hpp" ${CONCISE_ARGS_FILE} SHOW_PROGRESS)
endif()


add_executable(frsm-lcm frsm-lcm.cpp)
pods_use_pkg_config_packages(frsm-lcm frsm lcm bot2-core bot2-lcmgl-client)
pods_install_executables(frsm-lcm)




