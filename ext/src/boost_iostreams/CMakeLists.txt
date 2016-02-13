cmake_minimum_required(VERSION 2.8)

project(boost_iostreams CXX)

include_directories(${EXT_DIR}/include)

file(GLOB boost_iostreams_SRC
    "*.cpp"
)

add_definitions(-DNO_BZIP2=1)
add_library(boost_iostreams STATIC
            ${boost_iostreams_SRC})
