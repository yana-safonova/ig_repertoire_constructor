cmake_minimum_required(VERSION 2.8)

project(boost_filesystem CXX)

include_directories(${EXT_DIR}/include)

file(GLOB boost_filesystem_SRC
    "*.cpp"
)


add_definitions(-DBOOST_NOEXCEPT_OR_NOTHROW=noexcept)
add_library(boost_filesystem STATIC
            ${boost_filesystem_SRC})
