project(algorithms CXX)

INCLUDE_DIRECTORIES(SYSTEM ${SEQAN_INCLUDE_DIR})
include_directories(${IGREC_MAIN_INCLUDE_DIR})
include_directories(${CORE_DIR})
include_directories(${VDJ_UTILS_DIR})

add_library(algorithms STATIC
        block_alignment/block_alignment_primitives.cpp
        block_alignment/block_alignment_utils.cpp
        block_alignment/pairwise_block_alignment.cpp
        block_alignment/pairwise_block_aligner.cpp
        )

target_link_libraries(algorithms core vdj_utils)
