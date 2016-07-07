############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cdr_labeler CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})
include_directories(${CORE_DIR})
include_directories(${VDJ_UTILS_DIR})
include_directories(${ALGORITHMS_DIR})
include_directories(${VJ_FINDER_DIR})
include_directories(${GRAPH_UTILS})

add_executable(cdr_labeler
        cdr_config.cpp
        ../vj_finder/vj_finder_config.cpp
        ../vj_finder/germline_db_generator.cpp
        ../vj_finder/vj_alignment_structs.cpp
        ../vj_finder/vj_query_aligner.cpp
        ../vj_finder/vj_hits_filter.cpp
        ../vj_finder/vj_alignment_info.cpp
        ../vj_finder/vj_query_fix_fill_crop.cpp
        ../vj_finder/vj_query_processing.cpp
        ../vj_finder/vj_parallel_processor.cpp
        germline_db_labeler.cpp
        immunoglobulin_cdr_labeling/single_loop_labelers/single_loop_labeler.cpp
        immunoglobulin_cdr_labeling/single_loop_labelers/hcdr1_labeler.cpp
        immunoglobulin_cdr_labeling/single_loop_labelers/hcdr2_labeler.cpp
        immunoglobulin_cdr_labeling/single_loop_labelers/hcdr3_v_labeler.cpp
        immunoglobulin_cdr_labeling/single_loop_labelers/hcdr3_j_labeler.cpp
        immunoglobulin_cdr_labeling/single_loop_labelers/single_loop_helper.cpp
        immunoglobulin_cdr_labeling/immune_gene_labeler.cpp
        immunoglobulin_cdr_labeling/annotated_gene_labeler.cpp
        immunoglobulin_cdr_labeling/immune_gene_labeling_helper.cpp
        germline_db_labeling.cpp
        read_labeler.cpp
        compressed_cdr_set.cpp
        cdr_output.cpp
        diversity_analyser.cpp
        #cdr_annotator.cpp
        cdr_launch.cpp
        main.cpp)

target_link_libraries(cdr_labeler core algorithms vdj_utils graph_utils input boost_program_options ${COMMON_LIBRARIES})

