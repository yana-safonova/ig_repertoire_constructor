//
// Created by Andrew Bzikadze on 5/25/16.
//

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <logger/log_writers.hpp>

#include "../shm_kmer_model/evolutionary_edge_alignment/evolutionary_edge_alignment.hpp"
#include "../shm_kmer_model/alignment_checker/no_gaps_alignment_checker.hpp"
#include "../shm_kmer_model/shm_config.hpp"
#include "../shm_kmer_model/alignment_cropper/upto_last_reliable_kmer_alignment_cropper.hpp"
#include "../shm_kmer_model/germline_alignment_reader/alignment_reader.hpp"
#include "../shm_kmer_model/mutation_strategies/trivial_strategy.hpp"
#include "../shm_kmer_model/mutation_strategies/no_k_neighbours.hpp"

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

class AlignmentCheckerTest : public ::testing::Test {
public:
    void SetUp() { create_console_logger(); }
};

TEST_F(AlignmentCheckerTest, CheckingIsCorrect) {
    shm_config shm_config_ach_nogaps;
    std::string config_ach_nogaps = "test_dataset/shm_kmer_model/alignment_checker_nogaps.config.info";
    load(shm_config_ach_nogaps, config_ach_nogaps);
    NoGapsAlignmentChecker alch(shm_config_ach_nogaps.achp);

    ns_gene_alignment::EvolutionaryEdgeAlignment alignment1("AAA", "CCC", "id");
    ASSERT_TRUE(alch.check(alignment1));
    ns_gene_alignment::EvolutionaryEdgeAlignment alignment2("AAA--", "CCCAA", "id");
    ASSERT_FALSE(alch.check(alignment2));
    ns_gene_alignment::EvolutionaryEdgeAlignment alignment3("AAA--", "-CCAA", "id");
    ASSERT_FALSE(alch.check(alignment3));
}

class AlignmentCropperTest: public ::testing::Test {
public:
    void SetUp() { create_console_logger(); }
};

TEST_F(AlignmentCropperTest, CheckingIsCorrect) {
    shm_config shm_config_acr_last_rel_kmer;
    std::string config_acr_last_rel_kmer = "test_dataset/shm_kmer_model/alignment_cropper_last_rel_kmer.config.info";
    load(shm_config_acr_last_rel_kmer, config_acr_last_rel_kmer);
    UptoLastReliableKmerAlignmentCropper alcr(shm_config_acr_last_rel_kmer.acrp.rkmp);

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("AAAAAAAAAAAAAAA",
                                                           "AAAAAAAAAAAAAAA", "id");
        alcr.crop(alignment);
        ASSERT_EQ(alignment.son(),     "AAAAAAAAAAAAAAA");
        ASSERT_EQ(alignment.parent(), "AAAAAAAAAAAAAAA");
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("ABAAAAAAAAACAAA",
                                                           "ABAAAAAAAAACAAA", "id");
        alcr.crop(alignment);
        ASSERT_EQ(alignment.son(),     "ABAAAAAAAAACAAA");
        ASSERT_EQ(alignment.parent(), "ABAAAAAAAAACAAA");
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("BBAAAAAAAAACAAA",
                                                           "ABAAAAAAAAACAAA", "id");
        alcr.crop(alignment);
        ASSERT_EQ(alignment.son(),     "BAAAAAAAAACAAA");
        ASSERT_EQ(alignment.parent(), "BAAAAAAAAACAAA");
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("ABAAAAAAAAACAAA",
                                                           "ABAAAAAAAAACAAB", "id");
        alcr.crop(alignment);
        ASSERT_EQ(alignment.son(),     "ABAAAAAAAAACAA");
        ASSERT_EQ(alignment.parent(), "ABAAAAAAAAACAA");
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("ABAAAAAAAAACAAA",
                                                           "CBAAAAAAAAACAAB", "id");
        alcr.crop(alignment);
        ASSERT_EQ(alignment.son(),     "BAAAAAAAAACAA");
        ASSERT_EQ(alignment.parent(), "BAAAAAAAAACAA");
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("ABAAACAAAAACAAA",
                                                           "CBAAADAAAAACAAB", "id");
        alcr.crop(alignment);
        ASSERT_EQ(alignment.son(),     "AAAAACAA");
        ASSERT_EQ(alignment.parent(), "AAAAACAA");
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("ABAAACAAAAACAAA",
                                                           "CBAAADAAAAACAAB", "id");
        alcr.crop(alignment);
        ASSERT_EQ(alignment.son(),     "AAAAACAA");
        ASSERT_EQ(alignment.parent(), "AAAAACAA");
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("CAGGTGCAGCTGFTGC",
                                                           "CAGCCGCAGCTGGTGC", "id");
        alcr.crop(alignment);
        ASSERT_EQ(alignment.son(),     "GCAGCTG");
        ASSERT_EQ(alignment.parent(), "GCAGCTG");
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("CCGFTGCAGCTGFTTC",
                                                           "CAGhCGCAGCTGGTGC", "id");
        alcr.crop(alignment);
        ASSERT_EQ(alignment.son(),     "GCAGCTG");
        ASSERT_EQ(alignment.parent(), "GCAGCTG");
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("AAAAAA",
                                                           "AAAAAC", "id");
        alcr.crop(alignment);
        ASSERT_EQ(alignment.son(),     "AAAAA");
        ASSERT_EQ(alignment.parent(), "AAAAA");
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("CAGACGCAGCTGGTGF",
                                                           "AAGACGDAGCTGGTGC", "id");
        alcr.crop(alignment);
        ASSERT_EQ(alignment.son(),     "AGACGCAGCTGGTG");
        ASSERT_EQ(alignment.parent(), "AGACGDAGCTGGTG");
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("CAGGTGCAGZZGGTGCAGTCTGGGAGTGAACTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGCCGTCTATTACTGCACGAGACA",
                                                           "CAGGZGCAGCTGGTZCAGZCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCGGTCGTGTATTACTGTGZGAGAGA", "id");
        alcr.crop(alignment);
        ASSERT_EQ(alignment.son(),     "CTGGGAGTGAACTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGCCGTCTATTACTG");
        ASSERT_EQ(alignment.parent(), "CTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCGGTCGTGTATTACTG");
    }
}

class AlignmentReaderTest: public ::testing::Test {
public:
    void SetUp() { create_console_logger(); }
};


TEST_F(AlignmentReaderTest, CheckingIsCorrect) {
    shm_config shm_config_ar;
    std::string config_ar = "test_dataset/shm_kmer_model/germline_alignment_reader.config.info";

    load(shm_config_ar, config_ar);

    ns_alignment_reader::AlignmentReader alignment_reader(shm_config_ar.io.input.v_alignments,
                                                          shm_config_ar.io.input.cdr_details,
                                                          shm_config_ar.achp,
                                                          shm_config_ar.acrp);

    auto alignments = alignment_reader.read_alignments();
    ASSERT_EQ(alignments.size(), 4);
    ASSERT_EQ(alignments[0].son(), "CAGGTGCAGCTGGTGCAGTCTGGGAGTGAACTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGCCGTCTATTACTGCACGAGAGA");
    ASSERT_EQ(alignments[0].parent(), "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCGGTCGTGTATTACTGTGCGAGAGA");
    ASSERT_EQ(alignments[0].gene_id(), "cluster___8___size___1_HM855674|IGHV1-2*05|Homo sapiens|F|V-REGION|24..319|296 nt|1| | | | |296+0=296| | |");

    ASSERT_EQ(alignments[1].son(),     "CAGCTGCAGCTGCAGGAGTCGGGCCCAGGGCTGGTGAAGCCTTCGGGTCCAAGAACCAGTTCTCCCTAGAGGTGACCTCGGTGACCGCCGCAGACACGGCTGTGTATTAATGTGCGAG");
    ASSERT_EQ(alignments[1].parent(), "CAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCAGACACGGCTGTGTATTACTGTGCGAG");
    ASSERT_EQ(alignments[1].gene_id(),  "cluster___0___size___1_AB019439|IGHV4-39*01|Homo sapiens|F|V-REGION|11626..11924|299 nt|1| | | | |299+0=299| | |");

}

class MutationStrategiesTest: public ::testing::Test {
public:
    void SetUp() { create_console_logger(); }
};

TEST_F(MutationStrategiesTest, CheckNoKNeighbour) {
    shm_config shm_config_ms;
    std::string config_ms_nkn = "test_dataset/shm_kmer_model/mutation_strategy_nkn.config.info";
    load(shm_config_ms, config_ms_nkn);

    NoKNeighboursMutationStrategy ms_nkn(shm_config_ms.mfp);

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("AACCGGTTAA",
                                                           "AATCGGAAAA", "id");
        auto rel_pos = ms_nkn.calculate_relevant_positions(alignment);
        // ASSERT_THAT(v, ElementsAre(2, 5, 6, 7, 8, 9));
        ASSERT_EQ(rel_pos.size(), 1);
        ASSERT_EQ(rel_pos[0], 2);
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("AACCGGAAAA",
                                                           "AATCGGAAAA", "id");
        auto rel_pos = ms_nkn.calculate_relevant_positions(alignment);
        // ASSERT_THAT(v, ElementsAre(2, 5, 6, 7, 8, 9));
        ASSERT_EQ(rel_pos.size(), 4);
        ASSERT_EQ(rel_pos[0], 2);
        for (size_t i = 1; i < 4; ++i)
            ASSERT_EQ(rel_pos[i], i + 4);
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("TCTCCAACGTTTTCTG",
                                                           "CCTCCATCGGTTTCTG", "id");
        auto rel_pos = ms_nkn.calculate_relevant_positions(alignment);
        EXPECT_THAT(rel_pos, testing::ElementsAre(3, 6, 9, 12, 13));
    }

    {
        ns_gene_alignment::EvolutionaryEdgeAlignment alignment("TTTCCAACGTTTTCTGTGCACGAGGA",
                                                           "CCTCCATCGGTTTCTGTGCATGACGA", "id");
        auto rel_pos = ms_nkn.calculate_relevant_positions(alignment);
        EXPECT_THAT(rel_pos, testing::ElementsAre(6, 9, 12, 13, 14, 15, 16, 17, 20, 23));
    }
}
