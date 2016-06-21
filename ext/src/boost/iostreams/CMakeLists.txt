cmake_minimum_required(VERSION 2.8)

project(boost_iostreams CXX)

file(GLOB_RECURSE boost_iostreams_source_files "*.cpp")
add_library(boost_iostreams STATIC
            ${boost_iostreams_source_files})
