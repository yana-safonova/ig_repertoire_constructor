cmake_minimum_required(VERSION 2.8)

project(paired_read_merger CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${IGREC_MAIN_INCLUDE_DIR})
include_directories(${IGREC_MAIN_SRC_DIR}/ig_tools)

link_libraries(input yaml-cpp ${COMMON_LIBRARIES})
link_libraries(boost_iostreams)
add_executable(paired_read_merger main.cpp)
