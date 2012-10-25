
add_executable(frsm-test frsm-test.cpp)
pods_use_pkg_config_packages(frsm-test frsm)
target_link_libraries(frsm-test pthread)




