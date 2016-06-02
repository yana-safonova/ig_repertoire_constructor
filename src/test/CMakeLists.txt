project(tests)

include_directories(${CORE_DIR})
include_directories(${VDJ_UTILS_DIR})
include_directories(${ALGORITHMS_DIR})
include_directories(${VJ_FINDER_DIR})
include_directories(${CDR_LABELER_DIR})

link_libraries(gtest gmock_main gmock graph_utils vdj_utils algorithms input ${COMMON_LIBRARIES})

make_essential_test(test_sparse_graph test_sparse_graph.cpp)

# make_test(test_dsf test_dsf.cpp
#           ../dense_sgraph_finder/graph_decomposer/metis_permutation_constructor.cpp
#           ../dense_sgraph_finder/graph_decomposer/greedy_joining_decomposition_constructor.cpp
#           ../dense_sgraph_finder/graph_decomposer/simple_decomposition_constructor.cpp
#           ../dense_sgraph_finder/graph_decomposer/dense_subgraph_constructor.cpp)

add_executable(test_dsf test_dsf.cpp
               ../dense_sgraph_finder/graph_decomposer/metis_permutation_constructor.cpp
               ../dense_sgraph_finder/graph_decomposer/greedy_joining_decomposition_constructor.cpp
               ../dense_sgraph_finder/graph_decomposer/simple_decomposition_constructor.cpp
               ../dense_sgraph_finder/graph_decomposer/dense_subgraph_constructor.cpp)

add_executable(test_germline_database test_germline_database.cpp)

#add_executable(test_germline_cdr test_germline_cdr.cpp)

#add_executable(test_vjf_consistency test_vjf_consistency.cpp
#        ../vj_finder/vj_finder_config.cpp
#        ../vj_finder/germline_db_generator.cpp
#        ../vj_finder/vj_alignment_structs.cpp
#        ../vj_finder/vj_query_aligner.cpp
#        ../vj_finder/vj_hits_filter.cpp
#        ../vj_finder/vj_alignment_info.cpp
#        ../vj_finder/vj_query_fix_fill_crop.cpp
#        ../vj_finder/vj_query_processing.cpp)
