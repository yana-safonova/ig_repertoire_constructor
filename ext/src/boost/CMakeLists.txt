# -*- cmake -*-

cmake_minimum_required(VERSION 2.8)

add_definitions(-w)

include_directories(${EXT_DIR}/include)


add_subdirectory(filesystem)
add_subdirectory(iostreams)
add_subdirectory(program_options)
add_subdirectory(random)
add_subdirectory(regex)
add_subdirectory(system)
