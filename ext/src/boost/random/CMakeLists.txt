cmake_minimum_required(VERSION 2.8)

project(boost_random CXX)

file(GLOB_RECURSE boost_random_source_files "*.cpp")
add_library(boost_random STATIC
            ${boost_random_source_files})
