project(merged_reads_stats_calculator CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${IGREC_MAIN_INCLUDE_DIR})
include_directories(${IGREC_MAIN_SRC_DIR}/ig_tools)

link_libraries(boost_iostreams)
add_executable(merged_reads_stats_calculator
               main.cpp)
