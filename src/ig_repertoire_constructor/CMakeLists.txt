############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(ig_repertoire_constructor CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${CMAKE_SOURCE_DIR}/debruijn)

add_executable(ig_repertoire_constructor
               ../debruijn/kmer_coverage_model.cpp
               ../debruijn/config_struct.cpp
               ../debruijn/path_extend/pe_config_struct.cpp
               ../debruijn/stage.cpp
               ../debruijn/construction.cpp
               ../debruijn/gap_closer.cpp
               ../debruijn/simplification.cpp
               ../debruijn/mismatch_correction.cpp
               ../debruijn/pair_info_count.cpp
               ../debruijn/distance_estimation.cpp
               ../debruijn/repeat_resolving.cpp
               ../debruijn/genomic_info_filler.cpp
          read_archive.cpp
          ig_config.cpp          
          main.cpp)

target_link_libraries(ig_repertoire_constructor input mph_index BamTools yaml-cpp nlopt ${COMMON_LIBRARIES} ssw)

if (SPADES_STATIC_BUILD)
  set_target_properties(ig_repertoire_constructor PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS ig_repertoire_constructor
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/ig_repertoire_constructor"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")
