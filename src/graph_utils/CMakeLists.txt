cmake_minimum_required(VERSION 2.8)

file(GLOB graph_utils_SRC "*.cpp")

set(GRAPH_UTILS ${CMAKE_SOURCE_DIR}/graph_utils)

add_library(graph_utils STATIC ${graph_utils_SRC})
