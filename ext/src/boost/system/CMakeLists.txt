cmake_minimum_required(VERSION 2.8)

project(boost_system CXX)

file(GLOB_RECURSE boost_system_source_files "*.cpp")
add_library(boost_system STATIC
            ${boost_system_source_files})
