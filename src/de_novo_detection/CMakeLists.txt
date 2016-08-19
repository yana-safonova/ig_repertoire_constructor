############################################################################
# Copyright (c) 2011-2016 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(partial_order_graph CXX)

include_directories(SYSTEM ${SEQAN_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${PARTIAL_ORDER_GRAPH_DIR})

add_executable(de_novo_detection
    main.cpp
    )

target_link_libraries(de_novo_detection pog_library core boost_program_options ${COMMON_LIBRARIES})
