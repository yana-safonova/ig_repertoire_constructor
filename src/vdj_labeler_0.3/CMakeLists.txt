############################################################################
# Copyright (c) 2015-2015 Saint-Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(vdj_labeler_0.3 CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})
include_directories(${CORE_DIR})
include_directories(${ALGORITHMS_DIR})
include_directories(${VDJ_UTILS_DIR})
include_directories(${VJ_FINDER_DIR})
include_directories(SYSTEM ${SEQAN_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_executable(vdj_labeler_0.3
                vdj_config.cpp
                vdj_launch.cpp

                ../vj_finder/vj_finder_config.cpp
                ../vj_finder/germline_db_generator.cpp
                ../vj_finder/vj_alignment_structs.cpp
                ../vj_finder/vj_query_aligner.cpp
                ../vj_finder/vj_hits_filter.cpp
                ../vj_finder/vj_alignment_info.cpp
                ../vj_finder/vj_query_fix_fill_crop.cpp
                ../vj_finder/vj_query_processing.cpp
                ../vj_finder/vj_parallel_processor.cpp

                vdj_alignments/vdj_hits/single_immune_gene_segment_hits.cpp
                vdj_alignments/vdj_hits/vdj_hits.cpp
                vdj_alignments/vdj_hits/vdj_hits_storage.cpp
                vdj_alignments/vdj_hits/d_gene_segment_hit.cpp

                # vdj_alignments/aligners/simple_d_aligner.cpp
                # vdj_alignments/hits_calculator/alignment_quality_checkers/match_threshold_alignment_quality_checker.cpp
                # vdj_alignments/hits_calculator/d_hits_calculator/info_based_d_hits_calculator.cpp
                # vdj_alignments/hits_calculator/d_alignment_positions_calculator/custom_d_alignment_positions_calculator.cpp

                # model/recombination_model.cpp

                # recombination_calculator/hc_model_based_recombination_calculator.cpp
                # recombination_generation/gene_events_generators/shms_calculators/left_event_shms_calculator.cpp
                # recombination_generation/gene_events_generators/shms_calculators/right_event_shms_calculator.cpp
                # recombination_generation/gene_events_generators/shms_calculators/versatile_shms_calculator.cpp
                # recombination_generation/gene_events_generators/v_recombination_event_generator.cpp
                # recombination_generation/gene_events_generators/d_recombination_event_generator.cpp
                # recombination_generation/gene_events_generators/j_recombination_event_generator.cpp
                # recombination_generation/insertion_events_generators/versatile_insertion_event_generator.cpp
                # recombination_generation/custom_hc_recombination_generator.cpp
                # recombination_estimators/hc_recombination_estimator.cpp
                main.cpp
        )

target_link_libraries(vdj_labeler_0.3 core algorithms vdj_utils input yaml-cpp ${COMMON_LIBRARIES})
