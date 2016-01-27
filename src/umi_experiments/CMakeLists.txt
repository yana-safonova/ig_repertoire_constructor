project(umi_experiments CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

link_libraries(input ${COMMON_LIBRARIES} graph_utils)

add_executable(check_graph_symmetry check_graph_symmetry.cpp)
add_executable(umi_to_fastq umi_to_fastq.cpp)
add_executable(remove_stars remove_stars.cpp)
