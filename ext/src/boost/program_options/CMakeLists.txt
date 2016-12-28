cmake_minimum_required(VERSION 2.8)

project(boost_program_options CXX)

file(GLOB_RECURSE boost_program_options_source_files "*.cpp")
add_library(boost_program_options STATIC
            ${boost_program_options_source_files})
