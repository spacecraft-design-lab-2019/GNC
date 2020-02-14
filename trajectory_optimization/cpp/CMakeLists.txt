cmake_minimum_required(VERSION 3.13)
project(cpp)

set(CMAKE_CXX_STANDARD 14)

include_directories(../../eigen-git-mirror)

add_executable(iLQRsimple_test iLQRsimple.cpp PendulumTest.cpp utils.cpp)
add_executable(iLQRtest iLQR.cpp SatelliteTest.cpp utils.cpp iLQR.cpp)