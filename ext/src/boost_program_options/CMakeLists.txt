cmake_minimum_required(VERSION 2.8)

project(boost_program_options CXX)

include_directories(${EXT_DIR}/include)

file(GLOB boost_program_options_SRC
    "*.cpp"
)

add_library(boost_program_options STATIC
            ${boost_program_options_SRC})
