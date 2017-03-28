#include <gtest/gtest.h>

#include "../algorithms/block_alignment/block_alignment_utils.hpp"

using namespace algorithms;

TEST(block_chain_alignment_test, find_simple_gap_test) {
    EXPECT_EQ(1, find_simple_gap("A", "AC"));
    EXPECT_EQ(0, find_simple_gap("C", "AC"));
    EXPECT_EQ(1, find_simple_gap("AA", "ACA"));
    EXPECT_EQ(2, find_simple_gap("AC", "ACT"));
    EXPECT_EQ(0, find_simple_gap("CA", "ACA"));
    EXPECT_EQ(0, find_simple_gap("CAAAA", "ACAAAA"));
    EXPECT_EQ(3, find_simple_gap("ACAACA", "ACAXXXAXA"));
    EXPECT_EQ(0, find_simple_gap("", "ACAXXXAXA"));
    EXPECT_EQ(0, find_simple_gap("", "A"));
    EXPECT_EQ(0, find_simple_gap("", "GACTGGTTGG"));
}
