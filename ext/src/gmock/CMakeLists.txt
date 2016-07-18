cmake_minimum_required(VERSION 2.8)

project(gmock CXX)

include_directories(${EXT_DIR}/include)

add_library(gmock STATIC
            gmock-all.cc)

add_library(gmock_main STATIC
            gmock_main.cc)
