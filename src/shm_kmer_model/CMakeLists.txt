############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(shm_kmer_model CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})
include_directories(${GRAPH_UTILS})

add_executable(shm_kmer_model
               main.cpp shm_config.hpp shm_config.cpp shm_kmer_model_estimator.hpp shm_kmer_model_estimator.cpp alignment_reader.cpp gene_alignment.hpp alignment_checker/abstract_alignment_checker.hpp alignment_checker/no_gaps_alignment_checker.hpp alignment_checker/no_gaps_alignment_checker.cpp alignment_cropper/abstract_alignment_cropper.hpp alignment_cropper/upto_last_reliable_kmer_alignment_cropper.hpp alignment_cropper/upto_last_reliable_kmer_alignment_cropper.cpp mutations_strategies/abstract_mutation_strategy.hpp)

target_link_libraries(shm_kmer_model input yaml-cpp graph_utils ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(shm_kmer_model PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS shm_kmer_model
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/shm_kmer_model"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")
