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
#include "base_repertoire/metaroot/metaroot.hpp"
#include <seqan/seq_io.h>
#include "convert.hpp"

#include "base_repertoire/gene_chooser/uniform_gene_chooser.hpp"
#include "base_repertoire/nucleotides_remover/uniform_nucleotides_remover.hpp"
#include "base_repertoire/p_nucleotides_creator/uniform_nucleotides_creator.hpp"
#include "base_repertoire/n_nucleotides_inserter/uniform_n_nucleotides_inserter.hpp"
#include "base_repertoire/metaroot_creator/metaroot_creator.hpp"
#include "annotation_utils/cdr_labeling_primitives.hpp"
#include "base_repertoire/productivity_checker/productivity_checker.hpp"
#include "annotation_utils/aa_annotation/aa_calculator.hpp"

#include <chrono>
#include <random_generator.hpp>

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
        omp_set_num_threads(1);
        create_console_logger();
        std::string config_fname = "configs/ig_simulator/config.info";
        ig_simulator::load(config, config_fname);
        config.germline_params.loci = "IGH";

        germline_utils::GermlineDbGenerator db_generator(config.io_params.input_params.germline_input,
                                                         config.germline_params);
        v_db = db_generator.GenerateVariableDb();
        d_db = db_generator.GenerateDiversityDb();
        j_db = db_generator.GenerateJoinDb();
    }
};

TEST_F(IgSimulatorTest, PrepareGeneTest) {
    {
        seqan::Dna5String gene("GTACAACTGGAACG");
        AbstractMetaroot::PrepareGene(gene, 0, 1);
        ASSERT_EQ(core::seqan_string_to_string(gene), "GTACAACTGGAAC");
        AbstractMetaroot::PrepareGene(gene, 2, 1);
        ASSERT_EQ(core::seqan_string_to_string(gene), "ACAACTGGAA");
        AbstractMetaroot::PrepareGene(gene, 5, 3);
        ASSERT_EQ(core::seqan_string_to_string(gene), "TG");
        AbstractMetaroot::PrepareGene(gene, -2, -2);
        ASSERT_EQ(core::seqan_string_to_string(gene), "CATGCA");
        AbstractMetaroot::PrepareGene(gene, -3, 2);
        ASSERT_EQ(core::seqan_string_to_string(gene), "ATGCATG");
    }
}

TEST_F(IgSimulatorTest, VDJMetarootSequenceCorrect) {
    {
        std::string vd_ins("ACCGT");
        std::string dj_ins("TTTT");
        VDJMetaroot root(&v_db, &d_db, &j_db,
                         0, 0, 0,
                         annotation_utils::CDRLabeling(),
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
        VDJMetaroot root(&v_db, &d_db, &j_db,
                         0, 0, 0,
                         annotation_utils::CDRLabeling(),
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
        VDJMetaroot root(&v_db, &d_db, &j_db,
                         0, 0, 0,
                         annotation_utils::CDRLabeling(),
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
        VDJMetaroot root(&v_db, &d_db, &j_db,
                         0, 0, 0,
                         annotation_utils::CDRLabeling(),
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

TEST_F(IgSimulatorTest, VJMetarootSequenceCorrect) {
    {
        std::string vj_ins("ACCGT");
        VJMetaroot root(&v_db, &j_db,
                        0, 0,
                        annotation_utils::CDRLabeling(),
                        5, 3, vj_ins);
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

//TEST_F(IgSimulatorTest, MetarootCreaterSpeedTest) {
//    {
//        config.algorithm_params.germline_params.loci = "IGH";
//
//        germline_utils::GermlineDbGenerator db_generator(config.io_params.input_params.germline_input,
//                                                         config.algorithm_params.germline_params);
//        v_db = db_generator.GenerateVariableDb();
//        d_db = db_generator.GenerateDiversityDb();
//        j_db = db_generator.GenerateJoinDb();
//
//        VDJMetarootCreator metaroot_creator(config.simulation_params.base_repertoire_params.metaroot_simulation_params,
//                                            &v_db, &d_db, &j_db);
//
//        auto t1 = std::chrono::high_resolution_clock::now();
//        size_t N((int) 1e5);
//        for (size_t i = 0; i < N; ++i) {
//            auto root = metaroot_creator.Createroot();
//        }
//        auto t2 = std::chrono::high_resolution_clock::now();
//        std::chrono::duration<double, std::milli> fp = t2 - t1;
//        std::cout << "Simulation of " << N << " VDJ metaroots took " << fp.count() << "ms" << std::endl;
//    }
//
//    {
//        config.algorithm_params.germline_params.loci = "IGL";
//
//        germline_utils::GermlineDbGenerator db_generator(config.io_params.input_params.germline_input,
//                                                         config.algorithm_params.germline_params);
//        v_db = db_generator.GenerateVariableDb();
//        INFO("Generation of DB for join segments...");
//        j_db = db_generator.GenerateJoinDb();
//
//        VDJMetarootCreator metaroot_creator(config.simulation_params.base_repertoire_params.metaroot_simulation_params,
//                                            &v_db, &d_db, &j_db);
//
//        auto t1 = std::chrono::high_resolution_clock::now();
//        size_t N((int) 1e5);
//        for (size_t i = 0; i < N; ++i) {
//            metaroot_creator.Createroot()->Sequence();
//        }
//        auto t2 = std::chrono::high_resolution_clock::now();
//        std::chrono::duration<double, std::milli> fp = t2 - t1;
//        std::cout << "Simulation of " << N << " VJ metaroots took " << fp.count() << "ms" << std::endl;
//    }
//}

// TEST_F(IgSimulatorTest, MetarootCreaterCDRTest) {
//     {
//         config.algorithm_params.germline_params.loci = "IGH";
// 
//         germline_utils::GermlineDbGenerator db_generator(config.io_params.input_params.germline_input,
//                                                          config.algorithm_params.germline_params);
//         v_db = db_generator.GenerateVariableDb();
//         d_db = db_generator.GenerateDiversityDb();
//         j_db = db_generator.GenerateJoinDb();
// 
//         std::vector<germline_utils::CustomGeneDatabase> db;
//         db.emplace_back(std::move(v_db));
//         db.emplace_back(std::move(d_db));
//         db.emplace_back(std::move(j_db));
// 
//         VDJMetarootCreator metaroot_creator(config.simulation_params.base_repertoire_params.metaroot_simulation_params,
//                                             db);
// 
//         MTSingleton::SetSeed(7);
//         auto root = metaroot_creator.Createroot();
//         INFO(*root);
//         INFO(root->Sequence());
//         ASSERT_EQ(root->CDRLabeling().cdr1.start_pos, 75);
//         ASSERT_EQ(root->CDRLabeling().cdr1.end_pos, 98);
//         ASSERT_EQ(root->CDRLabeling().cdr2.start_pos, 150);
//         ASSERT_EQ(root->CDRLabeling().cdr2.end_pos, 173);
//         ASSERT_EQ(root->CDRLabeling().cdr3.start_pos, 288);
//         ASSERT_EQ(root->CDRLabeling().cdr3.end_pos, 366);
//     }
// }

TEST_F(IgSimulatorTest, ProductiveChecker) {
    {
        config.germline_params.loci = "IGH";

        germline_utils::GermlineDbGenerator db_generator(config.io_params.input_params.germline_input,
                                                         config.germline_params);
        v_db = db_generator.GenerateVariableDb();
        d_db = db_generator.GenerateDiversityDb();
        j_db = db_generator.GenerateJoinDb();

        std::vector<germline_utils::CustomGeneDatabase> db;
        db.emplace_back(std::move(v_db));
        db.emplace_back(std::move(d_db));
        db.emplace_back(std::move(j_db));

        ProductivityChecker productivity_checker(std::unique_ptr<annotation_utils::BaseAACalculator>
                                                (new annotation_utils::SimpleAACalculator));

        VDJMetarootCreator metaroot_creator(config.simulation_params.base_repertoire_params.metaroot_simulation_params,
                                            db);

        auto t1 = std::chrono::high_resolution_clock::now();
        size_t N((int) 1e4);
        size_t prod = 0;
        for (size_t i = 0; i < N; ++i) {
            auto root = metaroot_creator.Createroot();
            if (productivity_checker.IsProductive(root)) {
                prod++;
            }
        }
        std::cout << prod << " / " << N << std::endl;
        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> fp = t2 - t1;
        std::cout << "Simulation of " << N << " VDJ metaroots took " << fp.count() << "ms" << std::endl;
    }
}

} // End namespace ig_simulator
