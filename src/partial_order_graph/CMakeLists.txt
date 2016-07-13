############################################################################
# Copyright (c) 2011-2016 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(partial_order_graph CXX)

include_directories(SYSTEM ${SEQAN_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(partial_order_graph
    partial_order_graph.cpp
    pog_parameters.cpp
    pog_node.cpp
    main.cpp
    )

target_link_libraries(partial_order_graph core boost_program_options ${COMMON_LIBRARIES})
