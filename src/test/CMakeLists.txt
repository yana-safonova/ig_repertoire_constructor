project(tests)

include_directories(${CORE_DIR})
include_directories(${VDJ_UTILS_DIR})
include_directories(${ALGORITHMS_DIR})
include_directories(${VJ_FINDER_DIR})
include_directories(${CDR_LABELER_DIR})
include_directories(${SHM_DIR})
include_directories(${ANTEVOLO_DIR})
include_directories(${IG_SIMULATOR_DIR})
link_libraries(gtest gmock_main gmock graph_utils vdj_utils algorithms core input shm_kmer_matrix_estimator_library ${COMMON_LIBRARIES})

make_test(test_find_simple_gap test_find_simple_gap.cpp)
make_test(test_shm_kmer_model test_shm_kmer_model.cpp)
make_test(test_sparse_graph test_sparse_graph.cpp)
make_test(test_dsf test_dsf.cpp)

add_dependencies(test_dsf metis)
target_link_libraries(test_dsf dense_sgraph_finder_library)

make_test(test_sparse_graph test_sparse_graph.cpp)

make_test(test_antevolo test_antevolo.cpp
        ../antevolo/antevolo_config.cpp
        ../antevolo/shm_model_utils/shm_model.cpp
        ../antevolo/shm_model_utils/shm_model_edge_weight_calculator.cpp
        ../vj_finder/vj_finder_config.cpp
        ../vdj_utils/germline_utils/germline_db_generator.cpp
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
        ../cdr_labeler/read_labeler.cpp)

# make_test(test_dsf test_dsf.cpp)
# target_link_libraries(test_dsf dense_sgraph_finder_library)

# make_test(test_germline_database test_germline_database.cpp)

make_test(test_cdr_labeling test_cdr_labeler.cpp)
target_link_libraries(test_cdr_labeling cdr_labeler_library)

make_test(test_vj_finder test_vj_finder.cpp)
target_link_libraries(test_vj_finder vj_finder_library)

make_test(test_ig_simulator test_ig_simulator.cpp)
target_link_libraries(test_ig_simulator ig_simulator_library)
