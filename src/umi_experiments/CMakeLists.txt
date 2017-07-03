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
file(GLOB HEADER_FILES **/*.hpp)

#set(CMAKE_BUILD_TYPE "Debug")

add_executable(check_graph_symmetry tools/check_graph_symmetry.cpp ig_simulator_utils.cpp)
add_executable(umi_to_fastq tools/umi_to_fastq.cpp ig_simulator_utils.cpp umi_utils.cpp)
add_executable(simulate_tiny_dataset tools/simulate_tiny_dataset.cpp)
add_executable(umi_graph tools/umi_graph.cpp ig_simulator_utils.cpp umi_utils.cpp)
add_executable(analyze_intermed_clusters tools/analyze_intermed_clusters.cpp ig_simulator_utils.cpp)
add_executable(find_bad_cluster tools/find_bad_cluster.cpp ig_simulator_utils.cpp clusterer.cpp)
add_executable(report_pcr_error_rate tools/report_pcr_error_rate.cpp tools/error_analyzer.cpp tools/error_analyzer.hpp umi_utils.cpp utils/io.cpp)

add_executable(reads_by_umi_stats stats/reads_by_umi_stats.cpp ig_simulator_utils.cpp umi_utils.cpp stats/dist_distribution_stats.cpp)
add_executable(print_graph_decomposition_stats stats/print_graph_decomposition_stats.cpp ig_simulator_utils.cpp)
add_executable(pairwise_dist_stats stats/pairwise_dist_stats.cpp utils/io.cpp clusterer.cpp ig_simulator_utils.cpp)
add_executable(dists_inside_clusters stats/dists_inside_clusters.cpp utils/io.cpp clusterer.cpp ig_simulator_utils.cpp)
add_executable(umi_correction_stats stats/umi_correction_stats.cpp umi_utils.cpp utils/io.cpp)

add_executable(umi_naive naive/umi_naive.cpp ${HEADER_FILES} ig_simulator_utils.cpp umi_utils.cpp clusterer.cpp utils/io.cpp)
add_executable(report_umi_abundance report_umi_abundance.cpp utils/io.cpp)
add_executable(cluster_reads cluster_reads.cpp ${HEADER_FILES} ig_simulator_utils.cpp umi_utils.cpp clusterer.cpp utils/io.cpp ../fast_ig_tools/fast_ig_tools.cpp)
