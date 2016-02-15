project(tests)

link_libraries(gtest graph_utils input ${COMMON_LIBRARIES})

add_executable(test_sparse_graph test_sparse_graph.cpp)

add_executable(test_dsf test_dsf.cpp
        ../dense_sgraph_finder/graph_decomposer/metis_permutation_constructor.cpp
        ../dense_sgraph_finder/graph_decomposer/greedy_joining_decomposition_constructor.cpp
        ../dense_sgraph_finder/graph_decomposer/simple_decomposition_constructor.cpp
        ../dense_sgraph_finder/graph_decomposer/dense_subgraph_constructor.cpp)