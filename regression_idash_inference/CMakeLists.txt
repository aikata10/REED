project(HEregression)
cmake_minimum_required(VERSION 3.15)

find_package(SEAL 3.6 REQUIRED)

add_executable(regression_clean)

target_sources(regression_clean
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/main.cpp
)

target_link_libraries(regression_clean PRIVATE SEAL::seal)