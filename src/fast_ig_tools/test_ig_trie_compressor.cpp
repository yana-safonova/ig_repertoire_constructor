#include <gmock/gmock.h>

#include "ig_trie_compressor.hpp"

using namespace fast_ig_tools;
using namespace ::testing;

TEST(basic_tests, the_first_test) {
    std::vector<std::string> reads = {"AAA", "AAAA", "AAAC"};

    auto indices = compressed_reads_indices(reads);
    auto comp_reads = compressed_reads(reads.begin(), reads.end());

    EXPECT_THAT(indices.size(), reads.size());
    EXPECT_THAT(indices, ElementsAre(0, 0, 0));
    EXPECT_THAT(comp_reads.size(), 1);
    EXPECT_THAT(comp_reads, ElementsAre("AAA"));
}

TEST(basic_tests, empty_string_test) {
    std::vector<std::string> reads = {"AAA", "AAAA", "AAAC", "", "XXX", "sdadasdasd"};

    auto indices = compressed_reads_indices(reads);
    auto comp_reads = compressed_reads(reads);

    EXPECT_THAT(indices.size(), reads.size());
    EXPECT_THAT(indices, ElementsAre(3, 3, 3, 3, 3, 3));
    EXPECT_THAT(comp_reads.size(), 1);
    EXPECT_THAT(comp_reads, ElementsAre(""));
}

TEST(basic_tests, two_clusters) {
    std::vector<std::string> reads = {"AAA", "AAAA", "AAAC", "XX", "XXAA"};

    auto indices = compressed_reads_indices(reads);
    auto comp_reads = compressed_reads(reads);

    EXPECT_THAT(indices.size(), reads.size());
    EXPECT_THAT(indices, ElementsAre(0, 0, 0, 3, 3));
    EXPECT_THAT(comp_reads.size(), 2);
    EXPECT_THAT(comp_reads, ElementsAre("AAA", "XX"));
}
