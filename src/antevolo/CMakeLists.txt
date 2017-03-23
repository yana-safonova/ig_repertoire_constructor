############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(antevolo CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})
include_directories(${CORE_DIR})
include_directories(${VDJ_UTILS_DIR})
include_directories(${ALGORITHMS_DIR})
include_directories(${VJ_FINDER_DIR})
include_directories(${GRAPH_UTILS})
include_directories(${CDR_LABELER_DIR})
include_directories(${SHM_DIR})

add_executable(antevolo
        ../cdr_labeler/cdr_config.cpp
        ../vj_finder/vj_finder_config.cpp
        ../vj_finder/germline_db_generator.cpp
        ../vj_finder/vj_alignment_structs.cpp
        ../vj_finder/vj_query_aligner.cpp
        ../vj_finder/vj_hits_filter.cpp
        ../vj_finder/vj_alignment_info.cpp
        ../vj_finder/vj_query_fix_fill_crop.cpp
        ../vj_finder/vj_query_processing.cpp
        ../vj_finder/vj_parallel_processor.cpp
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
#        ../cdr_labeler/cdr_annotator.cpp
        ../cdr_labeler/compressed_cdr_set.cpp
        ../cdr_labeler/cdr_output.cpp
#        diversity_analyser.cpp
#        cdr_launch.cpp
#        ../fast_ig_tools/ig_trie_compressor.cpp
        clone_set_decomposers/vj_clone_set_decomposer.cpp
        vj_class_processors/clonally_related_candidates.cpp
        vj_class_processors/vj_class_processor.cpp
        vj_class_processors/edmonds_tarjan_DMST_calculator.cpp
        evolutionary_graph_utils/evolutionary_tree.cpp
        #evolutionary_graph_utils/evolutionary_graph.cpp
        #evolutionary_graph_utils/evolutionary_graph_constructor.cpp
        evolutionary_graph_utils/evolutionary_tree_splitter.cpp
        evolutionary_tree_storage.cpp
        evolutionary_tree_annotation/annotated_evolutionary_tree.cpp
        evolutionary_tree_annotation/tree_based_shm_annotator.cpp
        evolutionary_tree_annotation/evolutionary_annotated_shm.cpp
        evolutionary_tree_annotation/annotated_tree_storage.cpp
        evolutionary_stats_calculator.cpp
        antevolo_output_writer.cpp
        antevolo_config.cpp
        antevolo_processor.cpp
        antevolo_launch.cpp
        main.cpp
        cdr3_hamming_graph_connected_components_processors/base_cdr3_hg_cc_processor.cpp
        cdr3_hamming_graph_connected_components_processors/kruskal_cdr3_hg_cc_processor.cpp
        cdr3_hamming_graph_connected_components_processors/edmonds_cdr3_hg_cc_processor.cpp
        evolutionary_graph_utils/evolutionary_edge_constructor.cpp
        annotated_clone_by_read_constructor.cpp
        clone_set_with_fakes.hpp
        parent_read_reconstructor.cpp
        cdr3_hamming_graph_info.cpp
        related_clones_iterator.cpp
        shm_model_utils/shm_model.cpp
        posterior_distribution_calculator/posterior_distribution_calculator.cpp
        shm_model_utils/shm_model_edge_weight_calculator.cpp cluster_fillin_calculator.cpp cluster_fillin_calculator.hpp)

target_link_libraries(
        antevolo
        shm_kmer_matrix_estimator_library
        core
        algorithms
        vdj_utils
        graph_utils
        input
        boost_program_options
        ${COMMON_LIBRARIES})

