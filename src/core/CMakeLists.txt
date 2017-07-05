project(core CXX)

INCLUDE_DIRECTORIES(SYSTEM ${SEQAN_INCLUDE_DIR})
include_directories(${IGREC_MAIN_INCLUDE_DIR})

add_library(core STATIC
        read_archive.cpp
        aa_utils/amino_acid_motif.cpp
        aa_utils/aa_motif_finder.cpp
        convert.cpp
        )

target_link_libraries(core input)