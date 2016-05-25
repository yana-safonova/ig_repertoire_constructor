project(tests)

link_libraries(gtest gtest_main graph_utils input ${COMMON_LIBRARIES})

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
        ../shm_kmer_model/shm_config.hpp
        ../shm_kmer_model/shm_config.cpp
        ../shm_kmer_model/shm_kmer_model_estimator.hpp
        ../shm_kmer_model/shm_kmer_model_estimator.cpp
        ../shm_kmer_model/alignment_reader/alignment_reader.cpp
        ../shm_kmer_model/gene_alignment/gene_alignment.hpp
        ../shm_kmer_model/alignment_checker/abstract_alignment_checker.hpp
        ../shm_kmer_model/alignment_checker/no_gaps_alignment_checker.hpp
        ../shm_kmer_model/alignment_checker/no_gaps_alignment_checker.cpp
        ../shm_kmer_model/alignment_cropper/abstract_alignment_cropper.hpp
        ../shm_kmer_model/alignment_cropper/upto_last_reliable_kmer_alignment_cropper.hpp
        ../shm_kmer_model/alignment_cropper/upto_last_reliable_kmer_alignment_cropper.cpp
        ../shm_kmer_model/mutation_strategies/abstract_mutation_strategy.hpp
        ../shm_kmer_model/mutation_strategies/trivial_strategy.hpp
        ../shm_kmer_model/mutation_strategies/trivial_strategy.cpp
        ../shm_kmer_model/statistics_estimator/statistics_estimator.hpp
        ../shm_kmer_model/statistics_estimator/mutation_statistics.hpp
        ../shm_kmer_model/statistics_estimator/statistics_estimator.cpp
        ../shm_kmer_model/statistics_exporter/statistics_exporter.hpp
        ../shm_kmer_model/statistics_exporter/statistics_exporter.cpp
        ../shm_kmer_model/mutation_strategies/no_k_neighbours.hpp
        ../shm_kmer_model/mutation_strategies/no_k_neighbours.cpp
        ../shm_kmer_model/statistics_estimator/mutation_statistics.cpp)
