############################################################################
# Copyright (c) 2015-2015 Saint-Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(abpair CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})

find_path(SEQAN_INCLUDE_DIR seqan/seq_io.h
        HINTS ../../ext/include ${EXT_DIR}/include)

if(SEQAN_INCLUDE_DIR)
    message("seqan headers found in ${SEQAN_INCLUDE_DIR}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

INCLUDE_DIRECTORIES(SYSTEM ${SEQAN_INCLUDE_DIR})

add_executable(abpair
        cdr3_computation/simple_cdr3_calculator.cpp
        pairing_ig_data/pairing_primitives/droplet_barcode.cpp
        pairing_ig_data/pairing_primitives/ig_isotype.cpp
        pairing_ig_data/pairing_primitives/umi_isotype_sequences.cpp
        pairing_ig_data/pairing_fastq_utils.cpp
        pairing_ig_data/raw_pairing_data.cpp
        pairing_ig_data/db_factory.cpp
        pairing_ig_data/raw_pairing_data_storage.cpp
        pairing_ig_data/raw_pairing_data_stats_calculator.cpp
        abpair_config.cpp
        abpair_launch.cpp
        main.cpp)

target_link_libraries(abpair input ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
    set_target_properties(abpair PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS abpair
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/abpair"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")
