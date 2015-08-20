cmake_minimum_required(VERSION 2.8)

project(fast_ig_tools C CXX)

# Need newer cmake =(
# set_property(GLOBAL PROPERTY CXX_STANDARD 11)
# set(CMAKE_CXX_STANDARD 11)


include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
  add_definitions(-DSEQAN_HAS_ZLIB=1)
  link_libraries(${ZLIB_LIBRARIES})
else()
  add_definitions(-DSEQAN_HAS_ZLIB=0)
endif()

find_package(Boost REQUIRED program_options)

find_path(SEQAN_INCLUDE_DIR seqan/seq_io.h
          HINTS ../../ext/include ${EXT}/include)


find_library(TCMALLOC_LIBRARY NAMES tcmalloc)
if(TCMALLOC_LIBRARY)
  message("tcmalloc found: ${TCMALLOC_LIBRARY}")
endif()

find_library(JEMALLOC_LIBRARY NAMES jemalloc NO_DEFAULT_PATH PATHS ../../build/release/ext/jemalloc)
if(NOT JEMALLOC_LIBRARY)
  find_library(JEMALLOC_LIBRARY NAMES jemalloc HINT ../../build/release/ext/jemalloc)
endif()
if(JEMALLOC_LIBRARY)
  message("jemalloc found: ${JEMALLOC_LIBRARY}")
endif()

if(Boost_INCLUDE_DIR)
  message("Boost headers found in ${Boost_INCLUDE_DIR}")
endif()

if(SEQAN_INCLUDE_DIR)
  message("seqan headers found in ${SEQAN_INCLUDE_DIR}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

INCLUDE_DIRECTORIES(SYSTEM ${SEQAN_INCLUDE_DIR})
INCLUDE_DIRECTORIES(SYSTEM ${BOOST_INCLUDE_DIR})
link_libraries(pthread)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
find_library(BOOST_PROGRAM_OPTIONS_LIBRARY NAMES boost_program_options)
link_libraries(${BOOST_PROGRAM_OPTIONS_LIBRARY})

find_package(BZip2)
if(BZIP2_FOUND)
  add_definitions(-DSEQAN_HAS_BZIP2=1)
  link_libraries(${BZIP2_LIBRARIES})
else()
  add_definitions(-DSEQAN_HAS_BZIP2=0)
endif()

if(JEMALLOC_LIBRARY)
  message("Linking using jemalloc")
  link_libraries(${JEMALLOC_LIBRARY})
elseif(TCMALLOC_LIBRARY)
  message("Linking using tcmalloc")
  link_libraries(${TCMALLOC_LIBRARY})
else()
  message("Linking using ordinary malloc")
endif()

# Enable compilation in debug mode
# See http://stackoverflow.com/questions/19729036/boost-program-options-wont-work-with-glibcxx-debug
remove_definitions(-D_GLIBCXX_DEBUG)

add_executable(ig_trie_compressor ig_trie_compressor.cpp)
add_executable(ig_final_alignment ig_final_alignment.cpp)
add_executable(ig_kplus_vj_finder ig_kplus_vj_finder.cpp)
add_executable(ig_matcher ig_matcher.cpp)
add_executable(ig_swgraph_construct ig_swgraph_construct.cpp)
