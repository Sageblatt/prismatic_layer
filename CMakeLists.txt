cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11 -g -Og")

project(prismatic_layer)
add_executable(prism_layer main.cpp)

find_package(VTK REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

target_link_libraries(prism_layer Eigen3::Eigen)
target_link_libraries(prism_layer ${VTK_LIBRARIES})
