cmake_minimum_required(VERSION 3.15)
project(GNC)

set(CMAKE_CXX_STANDARD 14)

include_directories(orbit_propagation/cpp)
include_directories(util_funcs/cpp)

add_executable(GNC
        util_funcs/cpp/time_functions.cpp
        util_funcs/cpp/time_functions.h)
