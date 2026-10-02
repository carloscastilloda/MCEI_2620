cmake_minimum_required(VERSION 3.16)
project(Diff_Robot LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GSL REQUIRED)

add_executable(Diff_Robot main_diff_robot_CarlosCastillo.cpp)
target_link_libraries(Diff_Robot PRIVATE GSL::gsl GSL::gslcblas)
