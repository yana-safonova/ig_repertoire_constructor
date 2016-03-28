project(core CXX)

INCLUDE_DIRECTORIES(SYSTEM ${SEQAN_INCLUDE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})

add_library(core STATIC
        read_archive.cpp
        )

target_link_libraries(core input)