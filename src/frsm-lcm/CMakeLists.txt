pkg_check_modules(BOT_LCMGL bot2-lcmgl-client)
if(NOT BOT_LCMGL_FOUND)
    message("bot-lcmgl not found. Not building frsm-lcm executable")
    return()
endif()


#download the concise-args header
set(CONCISE_ARGS_FILE ${CMAKE_CURRENT_SOURCE_DIR}/ConciseArgs.hpp)
if(NOT EXISTS ${CONCISE_ARGS_FILE})
    file(DOWNLOAD "https://raw.githubusercontent.com/abachrach/concise-args/master/ConciseArgs.hpp" ${CONCISE_ARGS_FILE} SHOW_PROGRESS)
endif()


add_executable(frsm-lcm frsm-lcm.cpp)
pods_use_pkg_config_packages(frsm-lcm frsm lcm bot2-core bot2-lcmgl-client)
pods_install_executables(frsm-lcm)




