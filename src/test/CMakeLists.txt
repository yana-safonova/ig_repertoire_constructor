project(tests)

include_directories(${CORE_DIR})
include_directories(${VDJ_UTILS_DIR})
include_directories(${ALGORITHMS_DIR})
include_directories(${VJ_FINDER_DIR})
include_directories(${CDR_LABELER_DIR})
include_directories(${IG_SIMULATOR_DIR})

link_libraries(graph_utils vdj_utils algorithms core input ${COMMON_LIBRARIES})

make_test(test_find_simple_gap test_find_simple_gap.cpp)

make_test(test_sparse_graph test_sparse_graph.cpp)

make_test(test_dsf test_dsf.cpp)
add_dependencies(test_dsf metis)
target_link_libraries(test_dsf dense_sgraph_finder_library)

make_test(test_germline_database test_germline_database.cpp)

make_test(test_cdr_labeling test_cdr_labeler.cpp)
target_link_libraries(test_cdr_labeling cdr_labeler_library)

make_test(test_vj_finder test_vj_finder.cpp)
target_link_libraries(test_vj_finder vj_finder_library)

make_test(test_ig_simulator test_ig_simulator.cpp)
target_link_libraries(test_ig_simulator ig_simulator_library)
