project(umi_experiments CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

link_libraries(input ${COMMON_LIBRARIES} graph_utils boost_filesystem boost_system boost_program_options)

add_executable(check_graph_symmetry check_graph_symmetry.cpp)
add_executable(umi_to_fastq umi_to_fastq.cpp)
add_executable(print_graph_decomposition_stats print_graph_decomposition_stats.cpp graph_stats.hpp)
add_executable(umi_graph umi_graph.cpp)
