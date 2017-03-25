//
// Created by Andrew Bzikadze on 3/23/17.
//

#include <gtest/gtest.h>
#include <logger/log_writers.hpp>

#include <sstream>

#include <ig_simulator_config.hpp>
#include <germline_utils/germline_db_generator.hpp>
#include <germline_db_labeler.hpp>
#include <read_labeler.hpp>
#include "metaroot/metaroot.hpp"
#include <seqan/seq_io.h>
#include "convert.hpp"

#include "gene_chooser/uniform_gene_chooser.hpp"
#include "nucleotides_remover/uniform_nucleotides_remover.hpp"
#include "p_nucleotides_creator/uniform_nucleotides_creator.hpp"
#include "n_nucleotides_inserter/uniform_n_nucleotides_inserter.hpp"
#include "metaroot_creator/metaroot_creator.hpp"

#include <chrono>

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

ig_simulator::IgSimulatorConfig config;
germline_utils::CustomGeneDatabase v_db(germline_utils::SegmentType::VariableSegment);
germline_utils::CustomGeneDatabase d_db(germline_utils::SegmentType::DiversitySegment);
germline_utils::CustomGeneDatabase j_db(germline_utils::SegmentType::JoinSegment);

namespace ig_simulator {

class IgSimulatorTest: public ::testing::Test {
public:
    void SetUp() {
        create_console_logger();
        std::string config_fname = "configs/ig_simulator/config.info";
        ig_simulator::load(config, config_fname);
        config.algorithm_params.germline_params.loci = "IGH";

        germline_utils::GermlineDbGenerator db_generator(config.io_params.input_params.germline_input,
                                                         config.algorithm_params.germline_params);
        INFO("Generation of DB for variable segments...");
        v_db = db_generator.GenerateVariableDb();
        INFO("Generation of DB for diversity segments...");
        d_db = db_generator.GenerateDiversityDb();
        INFO("Generation of DB for join segments...");
        j_db = db_generator.GenerateJoinDb();
    }
};

TEST_F(IgSimulatorTest, PrepareGeneTest) {
    {
        seqan::Dna5String gene("GTACAACTGGAACG");
        AbstractMetaRoot::PrepareGene(gene, 0, 1);
        ASSERT_EQ(core::seqan_string_to_string(gene), "GTACAACTGGAAC");
        AbstractMetaRoot::PrepareGene(gene, 2, 1);
        ASSERT_EQ(core::seqan_string_to_string(gene), "ACAACTGGAA");
        AbstractMetaRoot::PrepareGene(gene, 5, 3);
        ASSERT_EQ(core::seqan_string_to_string(gene), "TG");
        AbstractMetaRoot::PrepareGene(gene, -2, -2);
        ASSERT_EQ(core::seqan_string_to_string(gene), "CATGCA");
        AbstractMetaRoot::PrepareGene(gene, -3, 2);
        ASSERT_EQ(core::seqan_string_to_string(gene), "ATGCATG");
    }
}

TEST_F(IgSimulatorTest, VDJMetaRootSequenceCorrect) {
    {
        std::string vd_ins("ACCGT");
        std::string dj_ins("TTTT");
        VDJMetaRoot root(&v_db, &d_db, &j_db,
                         0, 0, 0,
                         5, 1, 2, 3,
                         vd_ins, dj_ins);
        std::string correct_root_seq(
            std::string(
                "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTT"
                    "ACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACA"
                    "AACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGA"
                    "TCTGACGACACGGCCGTGTATTACTGTGCG") +
                vd_ins +
                "GTACAACTGGAACG" +
                dj_ins +
                "GAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG");
        ASSERT_EQ(correct_root_seq, core::dna5String_to_string(root.Sequence()));
    }

    {
        std::string vd_ins("ACCGT");
        std::string dj_ins("TTTT");
        VDJMetaRoot root(&v_db, &d_db, &j_db,
                         0, 0, 0,
                         -5, -2, -2, -3,
                         vd_ins, dj_ins);
        std::string correct_root_seq(
            std::string(
                "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTT"
                    "ACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACA"
                    "AACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGA"
                    "TCTGACGACACGGCCGTGTATTACTGTGCGAGAGA") +
                "TCTCT" +
                vd_ins +
                "CC" +
                "GGTACAACTGGAACGAC" +
                "GT" +
                dj_ins +
                "AGC" +
                "GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG");
        ASSERT_EQ(correct_root_seq, core::dna5String_to_string(root.Sequence()));
    }
    {
        std::string vd_ins("ACCGT");
        std::string dj_ins("TTTT");
        VDJMetaRoot root(&v_db, &d_db, &j_db,
                         0, 0, 0,
                         -5, 0, -3, -3,
                         vd_ins, dj_ins);
        std::string correct_root_seq(
            std::string(
                "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTT"
                    "ACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACA"
                    "AACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGA"
                    "TCTGACGACACGGCCGTGTATTACTGTGCGAGAGA") +
                "TCTCT" +
                vd_ins +
                "GGTACAACTGGAACGAC" +
                "GTC" +
                dj_ins +
                "AGC" +
                "GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG");
        ASSERT_EQ(correct_root_seq, core::dna5String_to_string(root.Sequence()));
    }

    {
        std::string vd_ins("ACCGT");
        std::string dj_ins("TTTT");
        VDJMetaRoot root(&v_db, &d_db, &j_db,
                         0, 0, 0,
                         0, 0, 0, 0,
                         vd_ins, dj_ins);
        std::string correct_root_seq(
            std::string(
                "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTT"
                    "ACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACA"
                    "AACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGA"
                    "TCTGACGACACGGCCGTGTATTACTGTGCGAGAGA") +
                vd_ins +
                "GGTACAACTGGAACGAC" +
                dj_ins +
                "GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG");
        ASSERT_EQ(correct_root_seq, core::dna5String_to_string(root.Sequence()));
    }
}

TEST_F(IgSimulatorTest, VJMetaRootSequenceCorrect) {
    {
        std::string vj_ins("ACCGT");
        VJMetaRoot root(&v_db, &j_db,
                        0, 0,
                        5, 3,
                        vj_ins);
        std::cout << root << std::endl;
        std::string correct_root_seq(
            std::string(
                "CAGGTTCAGCTGGTGCAGTCTGGAGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGTTACACCTTT"
                "ACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAGCGCTTACAATGGTAACACA"
                "AACTATGCACAGAAGCTCCAGGGCAGAGTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGA"
                "TCTGACGACACGGCCGTGTATTACTGTGCG") +
                vj_ins +
                "GAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG");
        ASSERT_EQ(correct_root_seq, core::dna5String_to_string(root.Sequence()));
    }
}

TEST_F(IgSimulatorTest, MetaRootCreater) {
    {
        AbstractVDJGeneChooserPtr gene_chooser(new UniformVDJGeneChooser(v_db, d_db, j_db));
        AbstractNucleotidesRemoverPtr nucl_remover(new UniformNucleotidesRemover());
        AbstractPNucleotidesCreatorPtr nucl_creator(new UniformPNucleotidesCreator());
        AbstractNNucleotidesInserterPtr nucl_inserter(new UniformNNucleotidesInserter());

        VDJMetarootCreator metaroot_creator(v_db, d_db, j_db,
                                            0.5, 0.5, 0.5, 0.5,
                                            std::move(gene_chooser),
                                            std::move(nucl_remover),
                                            std::move(nucl_creator),
                                            std::move(nucl_inserter));

        auto t1 = std::chrono::high_resolution_clock::now();
        size_t N((int) 1e6);
        for (size_t i = 0; i < N; ++i) {
            metaroot_creator.CreateRoot().Sequence();
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> fp = t2 - t1;
        std::cout << "Simulation of " << N << " VDJ metaroots took " << fp.count() << "ms" << std::endl;
    }

    {
        AbstractVDJGeneChooserPtr gene_chooser(new UniformVDJGeneChooser(v_db, j_db));
        AbstractNucleotidesRemoverPtr nucl_remover(new UniformNucleotidesRemover());
        AbstractPNucleotidesCreatorPtr nucl_creator(new UniformPNucleotidesCreator());
        AbstractNNucleotidesInserterPtr nucl_inserter(new UniformNNucleotidesInserter());

        VJMetarootCreator metaroot_creator(v_db, j_db,
                                           0.5, 0.5,
                                           std::move(gene_chooser),
                                           std::move(nucl_remover),
                                           std::move(nucl_creator),
                                           std::move(nucl_inserter));

        auto t1 = std::chrono::high_resolution_clock::now();
        size_t N((int) 1e6);
        for (size_t i = 0; i < N; ++i) {
            metaroot_creator.CreateRoot().Sequence();
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> fp = t2 - t1;
        std::cout << "Simulation of " << N << " VJ metaroots took " << fp.count() << "ms" << std::endl;
    }
}

} // End namespace ig_simulator