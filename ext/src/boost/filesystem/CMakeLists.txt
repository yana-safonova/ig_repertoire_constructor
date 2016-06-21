cmake_minimum_required(VERSION 2.8)

project(boost_filesystem CXX)

file(GLOB_RECURSE boost_filesystem_source_files "*.cpp")
add_library(boost_filesystem STATIC
            ${boost_filesystem_source_files})
