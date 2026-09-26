cmake_minimum_required( VERSION 3.20 )
project(integ_diff)
find_package( GSL )
add_executable( ${PROJECT_NAME} main_int_diff_CarlosCastillo.cpp )
