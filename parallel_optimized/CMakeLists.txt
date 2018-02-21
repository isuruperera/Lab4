cmake_minimum_required(VERSION 3.5)
project(Lab4)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

add_executable(MatrixMultiplyParallelOptimized main.cpp)
