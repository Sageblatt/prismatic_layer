cmake_minimum_required(VERSION 3.13)
project(prismatic_layer)

set(CMAKE_CXX_FLAGS "-std=c++11 -O2")

find_package(VTK REQUIRED)
find_package(Eigen3 REQUIRED)

set(SOURCES src)

add_library(prism_layer
        ${SOURCES}/Utils.cpp
        ${SOURCES}/Grid.cpp)

target_include_directories(prism_layer PRIVATE include)

target_link_libraries(prism_layer Eigen3::Eigen)
target_link_libraries(prism_layer ${VTK_LIBRARIES})

if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    add_executable(prism_layer_exec ${SOURCES}/main.cpp)
    target_include_directories(prism_layer_exec PRIVATE include)
    target_link_libraries(prism_layer_exec prism_layer)
    set_target_properties(prism_layer_exec PROPERTIES OUTPUT_NAME "prismatic_layer")
endif()