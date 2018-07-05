############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(antevolo CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${IGREC_MAIN_INCLUDE_DIR})
include_directories(${CORE_DIR})
include_directories(${VDJ_UTILS_DIR})
include_directories(${ALGORITHMS_DIR})
include_directories(${VJ_FINDER_DIR})
include_directories(${GRAPH_UTILS})
include_directories(${CDR_LABELER_DIR})
include_directories(${SHM_DIR})

add_library(antevolo_library STATIC
        clone_set_decomposers/vj_clone_set_decomposer.cpp
	    clone_set_decomposers/v_clone_set_decomposer.cpp
	    clone_set_decomposers/clone_set_decomposer.cpp
		vj_class_processors/vj_class_processor.cpp
        vj_class_processors/base_gene_class_processor.cpp
        #
        evolutionary_graph_utils/evolutionary_tree.cpp
        evolutionary_graph_utils/evolutionary_tree_splitter.cpp
        evolutionary_graph_utils/evolutionary_edge_constructor.cpp
        evolutionary_graph_utils/related_clones_iterator.cpp
        #
        evolutionary_tree_storage.cpp
        evolutionary_tree_annotation/annotated_evolutionary_tree.cpp
        evolutionary_tree_annotation/annotated_tree_storage.cpp
        #
        cdr3_hamming_graph_connected_components_processors/base_cdr3_hg_cc_processor.cpp
        cdr3_hamming_graph_connected_components_processors/kruskal_cdr3_hg_cc_processor.cpp
        cdr3_hamming_graph_connected_components_processors/edmonds_cdr3_hg_cc_processor.cpp
        cdr3_hamming_graph_connected_components_processors/edmonds_utils/edmonds_processor.cpp
        #
        shm_model_utils/shm_model.cpp
        posterior_distribution_calculator/posterior_distribution_calculator.cpp
        shm_model_utils/shm_model_edge_weight_calculator.cpp
		stats_calculation/cluster_fillin_calculator.cpp
		stats_calculation/cluster_fillin_calculator.hpp
        #
        shm_counting/tree_based_shm.cpp
        shm_counting/tree_shm_map.cpp
        shm_counting/tree_shm_calculator.cpp
        shm_counting/tree_based_shm_convertor.cpp
        #
        parallel_evolution/clonal_graph_writer.cpp
        parallel_evolution/clonal_graph.cpp
        parallel_evolution/clonal_graph_constructor.cpp
        parallel_evolution/aa_shm_map.cpp
        parallel_evolution/clonal_tree_refiner.cpp
        #
        annotated_clone_by_read_constructor.cpp
        clone_set_with_fakes.hpp
        parent_read_reconstructor.cpp
        cdr3_hamming_graph_component_info.cpp
        antevolo_output_writer.cpp
        antevolo_config.cpp
		reannotation/antevolo_reannotator.cpp
        antevolo_processor.cpp
        antevolo_launch.cpp
        evolutionary_graph_utils/one_child_fake_clones_filterer.cpp
		vj_class_processors/v_class_processor.cpp
		vj_class_processors/v_class_processor.hpp)

target_link_libraries(antevolo_library
        vj_finder_library
        cdr_labeler_library
        shm_kmer_matrix_estimator_library
        core
        algorithms
        vdj_utils
        graph_utils
        input
        boost_program_options
        ${COMMON_LIBRARIES})

add_executable(antevolo
        main.cpp)

target_link_libraries(antevolo
        antevolo_library
        )
