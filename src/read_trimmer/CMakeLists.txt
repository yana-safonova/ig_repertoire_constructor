############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(read_trimmer CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${IGREC_MAIN_INCLUDE_DIR})
include_directories(${CORE_DIR})

add_executable(read_trimmer
        main.cpp
        )

target_link_libraries(read_trimmer
        core
        build_info
        input
        boost_program_options
        ${COMMON_LIBRARIES}
        )
