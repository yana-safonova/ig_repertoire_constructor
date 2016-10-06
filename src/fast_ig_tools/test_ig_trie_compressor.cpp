#include <gmock/gmock.h>

#include "ig_trie_compressor.hpp"

using fast_ig_tools::Compressor;
using namespace ::testing;

TEST(basic_tests, the_first_test) {
    std::vector<std::string> reads = {"AAA", "AAAA", "AAAC"};

    auto indices = Compressor::compressed_reads_indices(reads, Compressor::Type::TrieCompressor);
    auto comp_reads = Compressor::compressed_reads(reads.begin(), reads.end(), Compressor::Type::TrieCompressor);

    EXPECT_THAT(indices.size(), reads.size());
    EXPECT_THAT(indices, ElementsAre(0, 0, 0));
    EXPECT_THAT(comp_reads, ElementsAre("AAA"));
}

TEST(basic_tests, empty_string_test) {
    std::vector<std::string> reads = {"AAA", "AAAA", "AAAC", "", "XXX", "sdadasdasd"};

    auto indices = Compressor::compressed_reads_indices(reads, Compressor::Type::TrieCompressor);
    auto comp_reads = Compressor::compressed_reads(reads, Compressor::Type::TrieCompressor);

    EXPECT_THAT(indices.size(), reads.size());
    EXPECT_THAT(indices, ElementsAre(3, 3, 3, 3, 3, 3));
    EXPECT_THAT(comp_reads, ElementsAre(""));
}

TEST(basic_tests, two_clusters) {
    std::vector<std::string> reads = {"AAA", "AAAA", "AAAC", "XX", "XXAA"};

    auto indices = Compressor::compressed_reads_indices(reads, Compressor::Type::TrieCompressor);
    auto comp_reads = Compressor::compressed_reads(reads, Compressor::Type::TrieCompressor);

    EXPECT_THAT(indices.size(), reads.size());
    EXPECT_THAT(indices, ElementsAre(0, 0, 0, 3, 3));
    EXPECT_THAT(comp_reads, ElementsAre("AAA", "XX"));
}

TEST(basic_tests, hashmap_test) {
    std::vector<std::string> reads = {"AAA", "AAAA", "AAA", "", "XXX", "sdadasdasd", "XXX", "X", "AAAA"};

    auto indices = Compressor::compressed_reads_indices(reads, Compressor::Type::HashCompressor);
    auto comp_reads = Compressor::compressed_reads(reads, Compressor::Type::HashCompressor);

    EXPECT_THAT(indices.size(), reads.size());
    EXPECT_THAT(indices, ElementsAre(0, 1, 0, 3, 4, 5, 4, 7, 1));
    EXPECT_THAT(comp_reads, ElementsAre("AAA", "AAAA", "", "XXX", "sdadasdasd", "X"));
}
