cmake_minimum_required(VERSION 2.8)

project(fast_ig_tools C CXX)

# Need newer cmake =(
# set_property(GLOBAL PROPERTY CXX_STANDARD 11)
# set(CMAKE_CXX_STANDARD 11)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})


# if(ZLIB_FOUND)
#   add_definitions(-DSEQAN_HAS_ZLIB=1)
#   link_libraries(${ZLIB_LIBRARIES})
# else()
#   add_definitions(-DSEQAN_HAS_ZLIB=0)
# endif()

add_definitions(-DSEQAN_HAS_ZLIB=1)

# This checking performed BEFORE compilation, so boost_po is not available at this moment
# find_package(Boost REQUIRED program_options)

find_path(SEQAN_INCLUDE_DIR seqan/seq_io.h
  HINTS ../../ext/include ${EXT_DIR}/include)

if(SEQAN_INCLUDE_DIR)
  message("seqan headers found in ${SEQAN_INCLUDE_DIR}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

INCLUDE_DIRECTORIES(SYSTEM ${SEQAN_INCLUDE_DIR})
# link_libraries(pthread)

link_libraries(boost_program_options)

# find_package(BZip2 QUIET)
# if(BZIP2_FOUND)
#   add_definitions(-DSEQAN_HAS_BZIP2=1)
#   link_libraries(${BZIP2_LIBRARIES})
# else()
#   add_definitions(-DSEQAN_HAS_BZIP2=0)
# endif()

# add_definitions(-DSEQAN_HAS_BZIP2=1)

# Enable compilation in debug mode
# See http://stackoverflow.com/questions/19729036/boost-program-options-wont-work-with-glibcxx-debug
remove_definitions(-D_GLIBCXX_DEBUG)

link_libraries(input yaml-cpp ${COMMON_LIBRARIES})

# TEMPORARY measure! Remove conversion warnings
remove_definitions(-Wconversion)

add_executable(ig_trie_compressor ig_trie_compressor.cpp)
target_link_libraries(ig_trie_compressor build_info)

add_executable(ig_kplus_vj_finder ig_kplus_vj_finder.cpp)
add_executable(ig_matcher ig_matcher.cpp)
add_executable(ig_swgraph_construct ig_swgraph_construct.cpp)
add_executable(ig_kmer_counter ig_kmer_counter.cpp)
add_executable(ig_hgc_complexity_estimator ig_hgc_complexity_estimator.cpp)
add_executable(ig_consensus_finder ig_consensus_finder.cpp)

# add_executable(test1 test_ig_kplus_vj_finder.cpp)
# target_link_libraries(test1 gtest_main gtest)
# add_test(NAME mytest1 COMMAND test1)
# add_test(NAME mytest2222 COMMAND test1)
# add_test(NAME mytest3 COMMAND test1)
#
# include(CTest)

# make_essential_test(test_ig_kplus_vj_finder test_ig_kplus_vj_finder.cpp)
# make_test(test1 test_ig_kplus_vj_finder.cpp)
# make_test(test2 test_ig_kplus_vj_finder.cpp)
# make_test(test3 test_ig_kplus_vj_finder.cpp)

# add_custom_target(meow COMMAND ${CMAKE_CTEST_COMMAND} -V)
# function(new_test test)
#   if (NOT TARGET ${test})
#     add_executable(${test} ${ARGN})
#     target_link_libraries(${test} gtest_main gtest)
#     add_test(NAME ${test} COMMAND ${test})
#   endif(NOT TARGET ${test})
#
#   add_dependencies(meow ${test})
#   set_target_properties(${test} PROPERTIES EXCLUDE_FROM_ALL 1)
# endfunction(make_test)



add_custom_command(TARGET ig_kplus_vj_finder POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${PROJECT_SOURCE_DIR}/germline"
        $<TARGET_FILE_DIR:ig_kplus_vj_finder>/germline)
