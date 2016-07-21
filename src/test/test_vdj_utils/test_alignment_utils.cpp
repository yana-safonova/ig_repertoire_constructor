//
// Created by Andrew Bzikadze on 7/8/16.
//

#include <gtest/gtest.h>
#include "alignment_utils/alignment_positions.hpp"
#include "alignment_utils/pairwise_alignment.hpp"
#include <seqan/align.h>
#include <iostream>

using namespace alignment_utils;
using namespace germline_utils;
using namespace seqan;
using namespace core;
using namespace std;

class AlignmentUtilsTest : public ::testing::Test { };

TEST_F(AlignmentUtilsTest, TestAlignmentPositions) {
    auto alignment_pos1 = AlignmentPositions();
    ASSERT_EQ(alignment_pos1.QueryAlignmentLength(), 1);
    ASSERT_EQ(alignment_pos1.SubjectAlignmentLength(), 1);

    auto alignment_pos2 = AlignmentPositions(make_pair<size_t, size_t>(1, 10),
                                             make_pair<size_t, size_t>(5, 15));
    ASSERT_EQ(alignment_pos2.QueryAlignmentLength(), 10);
    ASSERT_EQ(alignment_pos2.SubjectAlignmentLength(), 11);
}

TEST_F(AlignmentUtilsTest, TestImmuneGeneAlignmentPositions) {
    AlignmentPositions al_pos1(make_pair<size_t, size_t>(1, 10),
                               make_pair<size_t, size_t>(5, 15));
    ImmuneGene immune_gene;
    Read read;
    ImmuneGeneAlignmentPositions ig_alignment_pos1(al_pos1, immune_gene, read);
    EXPECT_EQ(ig_alignment_pos1.GeneStartPos(), 5);
    EXPECT_EQ(ig_alignment_pos1.GeneEndPos(), 15);
    EXPECT_EQ(ig_alignment_pos1.ReadStartPos(), 1);
    EXPECT_EQ(ig_alignment_pos1.ReadEndPos(), 10);
    EXPECT_EQ(ig_alignment_pos1.ReadAlignmentLength(), 10);
    EXPECT_EQ(ig_alignment_pos1.GeneAlignmentLength(), 11);
    EXPECT_FALSE(ig_alignment_pos1.IsEmpty());

    AlignmentPositions al_pos2(make_pair<size_t, size_t>(1, 0),
                               make_pair<size_t, size_t>(0, 5));
    ImmuneGeneAlignmentPositions ig_alignment_pos2(al_pos2, immune_gene, read);
    ASSERT_TRUE(ig_alignment_pos2.IsEmpty());
}

TEST_F(AlignmentUtilsTest, TestPairwiseAlignment) {
    Align<Dna5String, ArrayGaps> align;
    resize(rows(align), 2);
    string immune_gene_seq("AAA");
    string read_seq("GGG");
    ImmuneGene immune_gene(ImmuneGeneType(), "ImmuneGene", immune_gene_seq, 0);
    Read read("Read", read_seq, 0);
    {
        string seq1 = "ACCACCACAGT";
        string seq2 = "ACTTGCCACAGTA";
        assignSource(row(align, 0), seq1);
        assignSource(row(align, 1), seq2);
        double score = globalAlignment(align, Score<int, Simple>(0, -1, -1));

        ImmuneGeneReadAlignment igal(immune_gene, read, align, score);
        String<char, CStyle> query_str = igal.Query().seq;
        String<char, CStyle> subj_str = igal.Subject().seq();
        ASSERT_STREQ(query_str, read_seq.c_str());
        ASSERT_STREQ(subj_str, immune_gene_seq.c_str());
        EXPECT_EQ(igal.SubjectAlignmentLength(), 13);
        EXPECT_EQ(igal.QueryAlignmentLength(), 13);
        EXPECT_EQ(igal.AlignmentLength(), 13);
        EXPECT_EQ(igal.RealStartAlignmentPos(), 0);
        EXPECT_EQ(igal.RealEndAlignmentPos(), 11);
        EXPECT_EQ(igal.NumberMatches(), 9);
        EXPECT_EQ(igal.NumberMismatches(), 2);
        EXPECT_EQ(igal.NumberGaps(), 2);
        EXPECT_EQ(igal.NumberSHMs(), 4);
        EXPECT_DOUBLE_EQ(igal.Score(), score);
        EXPECT_DOUBLE_EQ(igal.NormalizedScore(), score / 13);

        EXPECT_EQ(igal.StartSubjectPosition(), 0);
        EXPECT_EQ(igal.EndSubjectPosition(), 11);
        EXPECT_EQ(igal.StartQueryPosition(), 0);
        EXPECT_EQ(igal.EndQueryPosition(), 12);
    }
    {
        string seq1 = "GGAAAGT";
        string seq2 = "AAAGT";
        assignSource(row(align, 0), seq1);
        assignSource(row(align, 1), seq2);
        double score = localAlignment(align, Score<int, Simple>(2, -1, -2));

        ImmuneGeneReadAlignment igal(immune_gene, read, align, score);
        EXPECT_EQ(igal.StartSubjectPosition(), 2);
        EXPECT_EQ(igal.StartQueryPosition(), 0);
        EXPECT_EQ(igal.EndSubjectPosition(), 6);
        EXPECT_EQ(igal.EndQueryPosition(), 4);
    }

    // TODO test on QueryPositionBySubjectPosition
}