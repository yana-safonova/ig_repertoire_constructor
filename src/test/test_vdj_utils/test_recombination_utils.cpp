//
// Created by Andrew Bzikadze on 7/9/16.
//

#include <gtest/gtest.h>
#include "recombination_utils/cleaved_gene.hpp"
#include "recombination_utils/nongenomic_insertion.hpp"
#include "recombination_utils/recombination.hpp"
#include "alignment_utils/pairwise_alignment.hpp"
#include "recombination_utils/recombination_storage.hpp"
#include "recombination_utils/insertion_event_storage.hpp"


using namespace recombination_utils;
using namespace alignment_utils;
using namespace germline_utils;
using namespace seqan;
using namespace core;
using namespace std;

class RecombinationUtilsTest : public ::testing::Test { };

TEST_F(RecombinationUtilsTest, TestCleavedGene) {
    Align<Dna5String, ArrayGaps> align;
    resize(rows(align), 2);
    ImmuneGene immune_gene;
    Read read;

    string seq1 = "ACCACCACAGT";
    string seq2 = "ACTTGCCACAGTA";
    assignSource(row(align, 0), seq1);
    assignSource(row(align, 1), seq2);
    double score = globalAlignment(align, Score<int, Simple>(0, -1, -1));
    cout << align;

    auto igal_ptr = make_shared<ImmuneGeneReadAlignment>(immune_gene, read, align, score);
    {
        CleavedIgGeneAlignment clal(igal_ptr, -1, 2, 1, 2);
        EXPECT_EQ(clal.GeneId(), size_t(-1));
        EXPECT_EQ(clal.SHMsNumber(), 1 + 2 + 4);

        auto clal2(clal);
        EXPECT_EQ(clal2.SHMsNumber(), 1 + 2 + 4);
        EXPECT_EQ(clal2.LeftCleavageLength(), -1);
        EXPECT_EQ(clal2.RightCleavageLength(), 2);
        EXPECT_TRUE(clal2.LeftEventIsPalindrome());
        EXPECT_FALSE(clal2.LeftEventIsCleavage());
        EXPECT_TRUE(clal2.RightEventIsCleavage());
        EXPECT_FALSE(clal2.RightEventIsPalindrome());
        EXPECT_EQ(clal2.GeneId(), size_t(-1));
        EXPECT_EQ(clal2.StartReadPosition(), size_t(-1));
        EXPECT_EQ(clal2.EndReadPosition(), 10);
    }
}

TEST_F(RecombinationUtilsTest, TestNongenomicInsertion) {
    NongenomicInsertion ni(3, 5);
    EXPECT_EQ(ni.length(), 3);
    EXPECT_TRUE(ni.Valid());

    NongenomicInsertion ni2(1, 0);
    EXPECT_EQ(ni2.length(), 0);
    EXPECT_TRUE(ni2.Valid());

    NongenomicInsertion ni3(2, 0);
    EXPECT_FALSE(ni3.Valid());
}

// TODO TestRecombination, TestStorages.
