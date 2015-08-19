cmake_minimum_required(VERSION 2.8)

project(fast_ig_tools CXX)

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

find_package(BZip2 REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Boost REQUIRED program_options)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

INCLUDE_DIRECTORIES(../../ext/include)

link_libraries(boost_program_options z bz2 pthread)

add_executable(ig_trie_compressor ig_trie_compressor.cpp)
add_executable(ig_final_alignment ig_final_alignment.cpp)
add_executable(ig_kplus_vj_finder ig_kplus_vj_finder.cpp)
add_executable(ig_matcher ig_matcher.cpp)
add_executable(ig_swgraph_construct ig_swgraph_construct.cpp)
