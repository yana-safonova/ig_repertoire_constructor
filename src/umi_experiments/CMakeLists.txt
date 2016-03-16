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

# needed for CLion to include them into the project
file(GLOB HEADER_FILES *.hpp)

set(CMAKE_BUILD_TYPE "Debug")

add_executable(check_graph_symmetry check_graph_symmetry.cpp utils.cpp)
add_executable(umi_to_fastq umi_to_fastq.cpp utils.cpp umi_utils.cpp)
add_executable(print_graph_decomposition_stats print_graph_decomposition_stats.cpp utils.cpp)
add_executable(umi_graph umi_graph.cpp utils.cpp umi_utils.cpp)
add_executable(group_by_umi group_by_umi.cpp utils.cpp umi_utils.cpp)
add_executable(reads_by_umi_stats reads_by_umi_stats.cpp utils.cpp umi_utils.cpp dist_distribution_stats.cpp)
add_executable(cluster_reads cluster_reads.cpp ${HEADER_FILES} utils.cpp umi_utils.cpp clusterer.cpp)
