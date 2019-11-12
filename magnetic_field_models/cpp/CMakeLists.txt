cmake_minimum_required(VERSION 3.13)
project(cpp)

set(CMAKE_CXX_STANDARD 14)

add_executable(cpp magnetic_field.cpp magnetic_field.h)

#include(FetchContent)
#FetchContent_Declare(
#        pybind11
#        GIT_REPOSITORY https://github.com/pybind/pybind11
#        GIT_TAG        v2.2.3
#)

#FetchContent_GetProperties(pybind11)
#if(NOT pybind11_POPULATED)
#    FetchContent_Populate(pybind11)
#    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
#endif()

#include_directories(orbit_propagation/cpp)
#include_directories(util_funcs/cpp)
#include_directories(eigen-git-mirror)

#pybind11_add_module(sun_utils_cpp util_funcs/cpp/sun_utils.cpp)