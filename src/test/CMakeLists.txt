project(tests)

include_directories(${CORE_DIR})
include_directories(${VDJ_UTILS_DIR})
include_directories(${ALGORITHMS_DIR})
include_directories(${VJ_FINDER_DIR})
include_directories(${CDR_LABELER_DIR})
include_directories(${SHM_DIR})

link_libraries(gtest gmock_main gmock graph_utils vdj_utils algorithms core input ${COMMON_LIBRARIES})

make_essential_test(test_sparse_graph test_sparse_graph.cpp)

# make_test(test_dsf test_dsf.cpp
#           ../dense_sgraph_finder/graph_decomposer/metis_permutation_constructor.cpp
#           ../dense_sgraph_finder/graph_decomposer/greedy_joining_decomposition_constructor.cpp
#           ../dense_sgraph_finder/graph_decomposer/simple_decomposition_constructor.cpp
#           ../dense_sgraph_finder/graph_decomposer/dense_subgraph_constructor.cpp)

# add_executable(test_dsf test_dsf.cpp
#                ../dense_sgraph_finder/graph_decomposer/metis_permutation_constructor.cpp
#                ../dense_sgraph_finder/graph_decomposer/greedy_joining_decomposition_constructor.cpp
#                ../dense_sgraph_finder/graph_decomposer/simple_decomposition_constructor.cpp
#                ../dense_sgraph_finder/graph_decomposer/dense_subgraph_constructor.cpp)

make_test(test_shm_kmer_model test_shm_kmer_model.cpp
        ../shm_kmer_matrix_estimator/shm_config.hpp
        ../shm_kmer_matrix_estimator/shm_config.cpp
        ../shm_kmer_matrix_estimator/shm_kmer_matrix_estimator.hpp
        ../shm_kmer_matrix_estimator/shm_kmer_matrix_estimator.cpp
        ../shm_kmer_matrix_estimator/germline_alignment_reader/alignment_reader.cpp
        ../shm_kmer_matrix_estimator/evolutionary_edge_alignment/evolutionary_edge_alignment.hpp
        ../shm_kmer_matrix_estimator/alignment_checker/abstract_alignment_checker.hpp
        ../shm_kmer_matrix_estimator/alignment_checker/no_gaps_alignment_checker.hpp
        ../shm_kmer_matrix_estimator/alignment_checker/no_gaps_alignment_checker.cpp
        ../shm_kmer_matrix_estimator/alignment_cropper/abstract_alignment_cropper.hpp
        ../shm_kmer_matrix_estimator/alignment_cropper/upto_last_reliable_kmer_alignment_cropper.hpp
        ../shm_kmer_matrix_estimator/alignment_cropper/upto_last_reliable_kmer_alignment_cropper.cpp
        ../shm_kmer_matrix_estimator/mutation_strategies/abstract_mutation_strategy.hpp
        ../shm_kmer_matrix_estimator/mutation_strategies/trivial_strategy.hpp
        ../shm_kmer_matrix_estimator/mutation_strategies/trivial_strategy.cpp
        ../shm_kmer_matrix_estimator/statistics_estimator/statistics_estimator.hpp
        ../shm_kmer_matrix_estimator/statistics_estimator/mutation_statistics.hpp
        ../shm_kmer_matrix_estimator/statistics_estimator/statistics_estimator.cpp
        ../shm_kmer_matrix_estimator/statistics_exporter/statistics_exporter.hpp
        ../shm_kmer_matrix_estimator/statistics_exporter/statistics_exporter.cpp
        ../shm_kmer_matrix_estimator/mutation_strategies/no_k_neighbours.hpp
        ../shm_kmer_matrix_estimator/mutation_strategies/no_k_neighbours.cpp
        ../shm_kmer_matrix_estimator/statistics_estimator/mutation_statistics.cpp)

add_executable(test_dsf test_dsf.cpp
               ../dense_sgraph_finder/graph_decomposer/metis_permutation_constructor.cpp
               ../dense_sgraph_finder/graph_decomposer/greedy_joining_decomposition_constructor.cpp
               ../dense_sgraph_finder/graph_decomposer/simple_decomposition_constructor.cpp
               ../dense_sgraph_finder/graph_decomposer/dense_subgraph_constructor.cpp)

add_executable(test_germline_database test_germline_database.cpp)

add_executable(test_cdr_labeling test_cdr_labeler.cpp
        ../vj_finder/vj_finder_config.cpp
        ../vj_finder/germline_db_generator.cpp
        ../vj_finder/vj_alignment_structs.cpp
        ../vj_finder/vj_query_aligner.cpp
        ../vj_finder/vj_hits_filter.cpp
        ../vj_finder/vj_alignment_info.cpp
        ../vj_finder/vj_query_fix_fill_crop.cpp
        ../vj_finder/vj_query_processing.cpp
        ../vj_finder/vj_parallel_processor.cpp
        ../cdr_labeler/cdr_config.cpp
        ../cdr_labeler/germline_db_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/single_loop_labelers/single_loop_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/single_loop_labelers/hcdr1_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/single_loop_labelers/hcdr2_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/single_loop_labelers/hcdr3_v_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/single_loop_labelers/hcdr3_j_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/single_loop_labelers/single_loop_helper.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/immune_gene_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/annotated_gene_labeler.cpp
        ../cdr_labeler/immunoglobulin_cdr_labeling/immune_gene_labeling_helper.cpp
        ../cdr_labeler/germline_db_labeling.cpp
        ../cdr_labeler/read_labeler.cpp
        )

add_executable(test_vj_finder test_vj_finder.cpp
        ../vj_finder/vj_finder_config.cpp
        ../vj_finder/germline_db_generator.cpp
        ../vj_finder/vj_alignment_structs.cpp
        ../vj_finder/vj_query_aligner.cpp
        ../vj_finder/vj_hits_filter.cpp
        ../vj_finder/vj_alignment_info.cpp
        ../vj_finder/vj_query_fix_fill_crop.cpp
        ../vj_finder/vj_query_processing.cpp
        ../vj_finder/vj_parallel_processor.cpp)
