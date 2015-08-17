INCLUDES=-isystem ./seqan-library-2.0.0/include/ -isystem ./gtest-1.7.0/fused-src/
CXXFLAGS=-std=c++11 -Wall -Wextra
CXXFLAGS+=-O3
CXXLIBS=-lpthread -lboost_program_options -lz -lbz2
CXXLIBS+=-lprofiler
CXXLIBS+=-ltcmalloc
CXXLIBS+=-fopenmp
# CXXLIBS+=-fopenmp=libiomp5
# CXXLIBS+=-ljemalloc
# CXX=clang++-3.6
CXX=g++

all: ig_swgraph_construct ig_trie_compressor ig_matcher ig_kplus_vj_finder ig_final_alignment

ig_swgraph_construct: ig_swgraph_construct.cpp ig_matcher.hpp
	${CXX} ${CXXFLAGS} ${INCLUDES} ig_swgraph_construct.cpp ${CXXLIBS} -o ig_swgraph_construct

ig_matcher: ig_matcher.cpp ig_matcher.hpp banded_half_smith_waterman.hpp
	${CXX} ${CXXFLAGS} ${INCLUDES} ig_matcher.cpp ${CXXLIBS} -o ig_matcher

ig_trie_compressor: ig_trie_compressor.cpp ig_trie_compressor.hpp
	${CXX} ${CXXFLAGS} ${INCLUDES} ig_trie_compressor.cpp ${CXXLIBS} -o ig_trie_compressor

ig_kplus_vj_finder: ig_kplus_vj_finder.cpp
	${CXX} ${CXXFLAGS} ${INCLUDES} ig_kplus_vj_finder.cpp ${CXXLIBS} -o ig_kplus_vj_finder

ig_final_alignment: ig_final_alignment.cpp ig_final_alignment.hpp
	${CXX} ${CXXFLAGS} ${INCLUDES} ig_final_alignment.cpp ${CXXLIBS} -o ig_final_alignment

gtest-all.o:
	${CXX} -c ${CXXFLAGS} ${INCLUDES} ./gtest-1.7.0/fused-src/gtest/gtest-all.cc -o gtest-all.o

gtest_main.o:
	${CXX} -c ${CXXFLAGS} ${INCLUDES} ./gtest-1.7.0/fused-src/gtest/gtest_main.cc -o gtest_main.o

clean:
	-rm *.o
	-rm ig_swgraph_construct
	-rm ig_kplus_vj_finder
	-rm ig_matcher
	-rm ig_trie_compressor
