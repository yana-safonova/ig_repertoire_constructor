project(umi_experiments CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

link_libraries(
        input
        ${COMMON_LIBRARIES}
        graph_utils
        boost_filesystem
        boost_system
        boost_program_options
)

set(SRC_FILES umi_utils.cpp dist_distribution_stats.cpp)

# needed for CLion to include them into the project
file(GLOB HEADER_FILES *.hpp)

add_executable(check_graph_symmetry check_graph_symmetry.cpp ${HEADER_FILES} ${SRC_FILES})
add_executable(umi_to_fastq umi_to_fastq.cpp ${SRC_FILES})
add_executable(print_graph_decomposition_stats print_graph_decomposition_stats.cpp ${SRC_FILES})
add_executable(umi_graph umi_graph.cpp ${SRC_FILES})
add_executable(group_by_umi group_by_umi.cpp ${SRC_FILES})
add_executable(reads_by_umi_stats reads_by_umi_stats.cpp ${SRC_FILES})
