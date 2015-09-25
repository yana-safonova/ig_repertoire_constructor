cmake_minimum_required(VERSION 2.8)

project(boost_filesystem CXX)

include_directories(${EXT_DIR}/include)

file(GLOB boost_filesystem_SRC
    "*.cpp"
)

add_library(boost_filesystem STATIC
            ${boost_filesystem_SRC})
