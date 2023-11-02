cmake_minimum_required(VERSION 3.2)
project(SubgraphMatching)

include_directories(
        ${PROJECT_SOURCE_DIR}/
)

set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -g -march=native -pthread -I/usr/local/opt/libomp/include")

add_subdirectory(graph)
add_subdirectory(matching)
add_subdirectory(utility)