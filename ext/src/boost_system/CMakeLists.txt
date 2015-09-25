cmake_minimum_required(VERSION 2.8)

project(boost_system CXX)

include_directories(${EXT_DIR}/include)

file(GLOB boost_system_SRC
    "*.cpp"
)

add_library(boost_system STATIC
            ${boost_system_SRC})
