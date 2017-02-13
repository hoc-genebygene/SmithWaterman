cmake_minimum_required(VERSION 3.0)
project(SmithWaterman CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Threads REQUIRED)

add_executable(SmithWaterman SmithWaterman.cpp)

