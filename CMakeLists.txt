cmake_minimum_required(VERSION 3.17)
project(PointPacking)

set(CMAKE_CXX_STANDARD 11)

add_executable(PointPacking
        main.cpp
        stat.cpp
        stat.h)