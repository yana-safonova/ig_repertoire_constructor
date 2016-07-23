//
// Created by Andrew Bzikadze on 7/13/16.
//

#include <gtest/gtest.h>
#include <iostream>
#include "alignment_utils/alignment_positions.hpp"
#include "alignment_utils/pairwise_alignment.hpp"
#include <seqan/align.h>
#include "recombination_generation/gene_events_generators/shms_calculators/left_event_shms_calculator.hpp"
#include "recombination_generation/gene_events_generators/shms_calculators/right_event_shms_calculator.hpp"

using namespace vdj_labeler;
using namespace alignment_utils;
using namespace germline_utils;
using namespace seqan;
using namespace core;
using namespace std;

class SHMCalculatorsTest: public ::testing::Test { };

TEST_F(SHMCalculatorsTest, left_event_shm_calculator) {
    Align<Dna5String, ArrayGaps> align;
    resize(rows(align), 2);
    germline_utils::ImmuneGeneType igt;
    LeftEventSHMsCalculator calculator;
    {
        ImmuneGene immune_gene(igt, "CoolImmuneGene", "CACCACCA", 0);
        Read read("EvenCoolerRead",               "GTTACACCACCA", 0);
        assignSource(row(align, 0), immune_gene.seq());
        assignSource(row(align, 1), read.seq);
        localAlignment(align, Score<int, Simple>(2, -1, -1));


        ImmuneGeneReadAlignment igal(immune_gene, read, align, 0);
        // cout << igal->Alignment();

        // Test zero
        EXPECT_EQ(calculator.ComputeNumberSHMsForLeftEvent(igal, 0), 0);

        // Test palindromes
        EXPECT_EQ(calculator.ComputeNumberSHMsForLeftEvent(igal, -1), 1);
        EXPECT_EQ(calculator.ComputeNumberSHMsForLeftEvent(igal, -2), 1);
        EXPECT_EQ(calculator.ComputeNumberSHMsForLeftEvent(igal, -3), 2);
        EXPECT_EQ(calculator.ComputeNumberSHMsForLeftEvent(igal, -3), 2);

        // Test cleavage
        for (int i = 0; i < static_cast<int>(igal.SubjectAlignmentLength()); ++i) {
            EXPECT_EQ(calculator.ComputeNumberSHMsForLeftEvent(igal, i), 0);
        }
    }
    {
        ImmuneGene immune_gene(igt, "CoolImmuneGene", "CACCACCA", 0);
        Read read("EvenCoolerRead",               "GTTACCCCACCA", 0);
        assignSource(row(align, 0), immune_gene.seq());
        assignSource(row(align, 1), read.seq);
        localAlignment(align, Score<int, Simple>(2, -1, -1));


        ImmuneGeneReadAlignment igal(immune_gene, read, align, 0);
        // cout << igal->Alignment();

        // Test zero
        EXPECT_EQ(calculator.ComputeNumberSHMsForLeftEvent(igal, 0), 0);

        // Test palindromes
        EXPECT_EQ(calculator.ComputeNumberSHMsForLeftEvent(igal, -1), 1);
        EXPECT_EQ(calculator.ComputeNumberSHMsForLeftEvent(igal, -2), 1);
        EXPECT_EQ(calculator.ComputeNumberSHMsForLeftEvent(igal, -3), 2);
        EXPECT_EQ(calculator.ComputeNumberSHMsForLeftEvent(igal, -3), 2);

        // Test cleavage
        EXPECT_EQ(calculator.ComputeNumberSHMsForLeftEvent(igal, 1), 0);
        for (int i = 2; i < static_cast<int>(igal.SubjectAlignmentLength()); ++i) {
            EXPECT_EQ(calculator.ComputeNumberSHMsForLeftEvent(igal, i), -1);
        }
    }
}

TEST_F(SHMCalculatorsTest, right_event_shm_calculator) {
    Align<Dna5String, ArrayGaps> align;
    resize(rows(align), 2);
    germline_utils::ImmuneGeneType igt;
    RightEventSHMsCalculator calculator;
    {
        ImmuneGene immune_gene(igt, "CoolImmuneGene", "ACCACCAC", 0);
        Read read("EvenCoolerRead",                   "ACCACCACATTG", 0);
        assignSource(row(align, 0), immune_gene.seq());
        assignSource(row(align, 1), read.seq);
        localAlignment(align, Score<int, Simple>(2, -1, -1));


        ImmuneGeneReadAlignment igal(immune_gene, read, align, 0);
        // cout << igal->Alignment();

        // Test zero
        EXPECT_EQ(calculator.ComputeNumberSHMsForRightEvent(igal, 0), 0);

        // Test palindromes
        EXPECT_EQ(calculator.ComputeNumberSHMsForRightEvent(igal, -1), 1);
        EXPECT_EQ(calculator.ComputeNumberSHMsForRightEvent(igal, -2), 1);
        EXPECT_EQ(calculator.ComputeNumberSHMsForRightEvent(igal, -3), 2);
        EXPECT_EQ(calculator.ComputeNumberSHMsForRightEvent(igal, -3), 2);

        // Test cleavage
        for (int i = 0; i < static_cast<int>(igal.SubjectAlignmentLength()); ++i) {
            EXPECT_EQ(calculator.ComputeNumberSHMsForRightEvent(igal, i), 0);
        }
    }
    {
        ImmuneGene immune_gene(igt, "CoolImmuneGene", "ACCACCAC", 0);
        Read read("EvenCoolerRead",                   "ACCACCCCATTG", 0);
        assignSource(row(align, 0), immune_gene.seq());
        assignSource(row(align, 1), read.seq);
        localAlignment(align, Score<double, Simple>(2., -0.5, -1.));


        ImmuneGeneReadAlignment igal(immune_gene, read, align, 0);
        // cout << igal->Alignment();

        // Test zero
        EXPECT_EQ(calculator.ComputeNumberSHMsForRightEvent(igal, 0), 0);

        // Test palindromes
        EXPECT_EQ(calculator.ComputeNumberSHMsForRightEvent(igal, -1), 1);
        EXPECT_EQ(calculator.ComputeNumberSHMsForRightEvent(igal, -2), 1);
        EXPECT_EQ(calculator.ComputeNumberSHMsForRightEvent(igal, -3), 2);
        EXPECT_EQ(calculator.ComputeNumberSHMsForRightEvent(igal, -3), 2);

        // Test cleavage
        EXPECT_EQ(calculator.ComputeNumberSHMsForRightEvent(igal, 1), 0);
        for (int i = 2; i < static_cast<int>(igal.SubjectAlignmentLength()); ++i) {
            EXPECT_EQ(calculator.ComputeNumberSHMsForRightEvent(igal, i), -1);
        }
    }
}
