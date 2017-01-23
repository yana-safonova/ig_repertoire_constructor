############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(shm_kmer_matrix_estimator CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})
include_directories(${GRAPH_UTILS})

add_executable(shm_kmer_matrix_estimator
        main.cpp
        shm_config.cpp
        shm_kmer_matrix_estimator.cpp
        germline_alignment_reader/alignment_reader.cpp
        alignment_checker/no_gaps_alignment_checker.cpp
        alignment_cropper/upto_last_reliable_kmer_alignment_cropper.cpp
        mutation_strategies/trivial_strategy.cpp
        mutation_strategies/no_k_neighbours.cpp
        statistics_estimator/statistics_estimator.cpp
        statistics_estimator/mutation_statistics.cpp
        statistics_exporter/statistics_exporter.cpp
        command_line_routines.hpp
        command_line_routines.cpp
        )

target_link_libraries(shm_kmer_matrix_estimator boost_program_options input yaml-cpp graph_utils ${COMMON_LIBRARIES})
