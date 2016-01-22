############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(dense_sgraph_finder CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})
include_directories(${GRAPH_UTILS})

add_executable(dense_sgraph_finder
               graph_decomposer/metis_permutation_constructor.cpp
               graph_decomposer/greedy_joining_decomposition_constructor.cpp
               graph_decomposer/simple_decomposition_constructor.cpp
               graph_decomposer/dense_subgraph_constructor.cpp
               graph_decomposer/decomposition_stats_calculator.cpp
               dsf_config.cpp
               launch.cpp
               main.cpp)

target_link_libraries(dense_sgraph_finder input yaml-cpp graph_utils ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(dense_sgraph_finder PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS dense_sgraph_finder
        DESTINATION bin
        COMPONENT runtime)
install(DIRECTORY "${SPADES_CFG_DIR}/dense_sgraph_finder"
        DESTINATION share/spades/configs
        FILES_MATCHING PATTERN "*.info.template")
