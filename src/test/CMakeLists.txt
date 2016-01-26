project(tests)

link_libraries(gtest graph_utils input ${COMMON_LIBRARIES})

add_executable(test_sparse_graph test_sparse_graph.cpp)