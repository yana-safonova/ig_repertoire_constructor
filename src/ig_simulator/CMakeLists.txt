cmake_minimum_required(VERSION 2.8.4)
project(ig_simulator)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})
include_directories(${GRAPH_UTILS})

add_executable(ig_simulator
        main.cpp productivity_checker.hpp)

if (SPADES_STATIC_BUILD)
    set_target_properties(shm_kmer_model PROPERTIES LINK_SEARCH_END_STATIC 1)
endif ()
