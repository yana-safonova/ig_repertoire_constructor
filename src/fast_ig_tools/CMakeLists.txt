cmake_minimum_required(VERSION 2.8)

project(fast_ig_tools CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${SPADES_MAIN_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Enable compilation in debug mode
# See http://stackoverflow.com/questions/19729036/boost-program-options-wont-work-with-glibcxx-debug
remove_definitions(-D_GLIBCXX_DEBUG)

link_libraries(input yaml-cpp ${COMMON_LIBRARIES})
link_libraries(boost_program_options)

add_executable(ig_trie_compressor ig_trie_compressor.cpp fast_ig_tools.cpp utils.cpp)
target_link_libraries(ig_trie_compressor build_info)
target_link_libraries(ig_trie_compressor boost_system)

add_executable(ig_matcher ig_matcher.cpp fast_ig_tools.cpp utils.cpp)

add_executable(ig_swgraph_construct ig_swgraph_construct.cpp fast_ig_tools.cpp utils.cpp)
target_link_libraries(ig_swgraph_construct build_info)

add_executable(ig_kmer_counter ig_kmer_counter.cpp fast_ig_tools.cpp utils.cpp)
add_executable(ig_hgc_complexity_estimator ig_hgc_complexity_estimator.cpp fast_ig_tools.cpp utils.cpp)
add_executable(ig_consensus_finder ig_consensus_finder.cpp fast_ig_tools.cpp utils.cpp)

add_executable(ig_component_splitter ig_component_splitter.cpp utils.cpp)
target_link_libraries(ig_component_splitter build_info)

make_test(test_ig_trie_compressor test_ig_trie_compressor.cpp)
