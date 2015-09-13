cmake_minimum_required(VERSION 2.8)

project(boost_regex CXX)

include_directories(${EXT_DIR}/include)

file(GLOB boost_regex_SRC
    "*.cpp"
)

add_library(boost_regex STATIC
            ${boost_regex_SRC})
