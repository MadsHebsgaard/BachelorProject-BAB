cmake_minimum_required(VERSION 3.22)
project(Bachelor)

set(CMAKE_CXX_STANDARD 23)

add_executable(Bachelor main.cpp load.h save.h ui.h calculations.h testing.h createData.h load_output.h)
