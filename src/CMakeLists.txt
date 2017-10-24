# -*- cmake -*-
############################################################################
# Copyright (c) 2011-2014 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################
project(src)
# Hard prerequisites
cmake_minimum_required(VERSION 2.8.8)

# export compiler arguments (need for autocompletion)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# allow more human readable "if then else" constructs
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

# modules to find additional libraries (like zlib, boost and etc)
# Add path for custom modules
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")
# Define various dirs

set(IGREC_MAIN_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(IGREC_MAIN_INCLUDE_DIR ${IGREC_MAIN_SRC_DIR}/include)
set(IGREC_BINARY_DIR ${CMAKE_BINARY_DIR})
set(IGREC_TOOLS_BINARY_DIR ${IGREC_BINARY_DIR}/bin)
set(IGREC_BUILT_INCLUDE_DIR ${IGREC_BINARY_DIR}/include)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${IGREC_TOOLS_BINARY_DIR})
set(EXT_DIR "${IGREC_MAIN_SRC_DIR}/../ext")
set(IGREC_CFG_DIR "${IGREC_MAIN_SRC_DIR}/../configs")

set(CORE_DIR "${IGREC_MAIN_SRC_DIR}/core")
set(ALGORITHMS_DIR "${IGREC_MAIN_SRC_DIR}/algorithms")
set(VDJ_UTILS_DIR "${IGREC_MAIN_SRC_DIR}/vdj_utils")
set(VJ_FINDER_DIR "${IGREC_MAIN_SRC_DIR}/vj_finder")
set(CDR_LABELER_DIR "${IGREC_MAIN_SRC_DIR}/cdr_labeler")
set(IG_SIMULATOR_DIR "${IGREC_MAIN_SRC_DIR}/ig_simulator")
set(SHM_DIR "${IGREC_MAIN_SRC_DIR}/shm_kmer_matrix_estimator")
set(ANTEVOLO_DIR "${IGREC_MAIN_SRC_DIR}/antevolo")

# Everything option-dependent
include(options)

# Sanity check our source directory to make sure that we are not trying to
# generate an in-tree build
if ((CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR) OR
    (CMAKE_SOURCE_DIR STREQUAL "${CMAKE_BINARY_DIR}/src"))
  message(FATAL_ERROR "In-source builds are not allowed.
Please create a directory and run cmake from there, passing the path
to this source directory as the last argument.
This process created the file `CMakeCache.txt' and the directory `CMakeFiles'.
Please delete them.")
endif()

# Handle commond deps
include(deps)

# Add flags if necessary
include(flags)

# Handle common includes
include(includes)

# Handle common libs
include(libs)

# Enable testing (check & check_essential targets)
include(check)

configure_file("${IGREC_MAIN_INCLUDE_DIR}/config.hpp.in"
               "${IGREC_BUILT_INCLUDE_DIR}/config.hpp")

# Build external dependencies (if any)
add_subdirectory("${EXT_DIR}/src" "${CMAKE_BINARY_DIR}/ext")
add_subdirectory("${EXT_DIR}/tools" "${CMAKE_BINARY_DIR}/ext_tools")

add_subdirectory(build_info)

# sub projects
add_subdirectory(io)
add_subdirectory(fast_ig_tools)
add_subdirectory(dense_sgraph_finder)
add_subdirectory(paired_read_merger)
add_subdirectory(ig_tools)
add_subdirectory(graph_utils)
add_subdirectory(shm_kmer_matrix_estimator)
add_subdirectory(core)
add_subdirectory(algorithms)
add_subdirectory(vdj_utils)
add_subdirectory(vj_finder)
add_subdirectory(cdr_labeler)
add_subdirectory(antevolo)
add_subdirectory(test)
add_subdirectory(umi_experiments)
add_subdirectory(pcr_simulator)
add_subdirectory(ig_simulator)
add_subdirectory(config)
