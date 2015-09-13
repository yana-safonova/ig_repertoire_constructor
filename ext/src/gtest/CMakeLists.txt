cmake_minimum_required(VERSION 2.8)

project(gtest CXX)

include_directories(${EXT_DIR}/include)

add_library(gtest STATIC
            gtest-all.cc)
add_library(gtest_main STATIC
            gtest_main.cc)
