############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(shm_kmer_matrix_estimator CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${IGREC_MAIN_INCLUDE_DIR})
include_directories(${GRAPH_UTILS})
include_directories(${VDJ_UTILS_DIR})
include_directories(${CORE_DIR})

add_library(shm_kmer_matrix_estimator_library STATIC
        germline_alignment_reader/alignment_reader.cpp
        alignment_checker/no_gaps_alignment_checker.cpp
        alignment_cropper/upto_last_reliable_kmer_alignment_cropper.cpp
        mutation_strategies/trivial_strategy.cpp
        mutation_strategies/no_k_neighbours.cpp
        kmer_utils/kmer_utils.cpp
        kmer_utils/kmer_indexed_vector.cpp
        shm_kmer_matrix_estimator/shm_kmer_matrix_estimator.cpp
        kmer_matrix_exporter/kmer_matrix_exporter.cpp
        shm_kmer_matrix_estimator_config.cpp
        )


add_executable(shm_kmer_matrix_estimator
        main.cpp
        shm_kmer_matrix_estimator_pipeline.cpp
        )

target_link_libraries(shm_kmer_matrix_estimator
                      shm_kmer_matrix_estimator_library
                      vdj_utils
                      core
                      boost_program_options
                      input
                      yaml-cpp
                      graph_utils
                      ${COMMON_LIBRARIES})
