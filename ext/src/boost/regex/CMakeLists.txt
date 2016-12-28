cmake_minimum_required(VERSION 2.8)

project(boost_regex CXX)

file(GLOB_RECURSE boost_regex_source_files "*.cpp")
add_library(boost_regex STATIC
            ${boost_regex_source_files})
