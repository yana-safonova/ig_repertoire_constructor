############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(vj_finder CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})
include_directories(${CORE_DIR})
include_directories(${VDJ_UTILS_DIR})
include_directories(${ALGORITHMS_DIR})

add_definitions(-DSEQAN_HAS_ZLIB=1)
add_definitions(-DSEQAN_HAS_BZIP2=1)

add_executable(vj_finder
        vj_finder_config.cpp
        command_line_routines.cpp
        germline_db_generator.cpp
        vj_alignment_structs.cpp
        vj_query_aligner.cpp
        vj_hits_filter.cpp
        vj_alignment_info.cpp
        vj_query_fix_fill_crop.cpp
        vj_query_processing.cpp
        vj_parallel_processor.cpp
        vjf_launch.cpp
        main.cpp
        )

target_link_libraries(vj_finder core algorithms vdj_utils build_info input boost_program_options ${COMMON_LIBRARIES})

