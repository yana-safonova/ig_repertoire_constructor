//
// Created by Yana Safonova on 23/03/16.
//

#include <gtest/gtest.h>
#include <logger/log_writers.hpp>

#include <germline_utils/germline_databases/immune_gene_database.hpp>
#include <germline_utils/germline_databases/chain_database.hpp>
#include <germline_utils/germline_databases/custom_gene_database.hpp>

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

class GermlineDBTest : public ::testing::Test {
public:
    void SetUp() {
        create_console_logger();
    }
};

// test creates DB for specific type of immune gene (IGHV), uploads sequences from file and
// check correctness of the created DB
TEST_F(GermlineDBTest, TestImmuneGeneDbCorrectness) {
    std::string human_ighv_file = "data/germline/human/IG/IGHV.fa";
    using namespace germline_utils;
    ImmuneGeneDatabase ighv_database(ImmuneGeneType(ChainType(ImmuneChainType::HeavyIgChain),
                                                    SegmentType::VariableSegment));
    INFO("Creating database for immune gene " << ighv_database.GeneType());
    ighv_database.AddGenesFromFile(human_ighv_file);
    ASSERT_EQ(ighv_database.size(), 315);
    INFO("DB contains " << ighv_database.size() << " records");
    std::stringstream chain_str;
    chain_str << ighv_database.Chain();
    ASSERT_EQ(chain_str.str(), "IGH");
    INFO("DB chain type: " << ighv_database.Chain());
    std::stringstream gene_str;
    gene_str << ighv_database.GeneType();
    ASSERT_EQ(gene_str.str(), "IGHV");
    INFO("DB gene type: " << ighv_database.GeneType());
    //ASSERT_EQ(ighv_database[1].name(), "IGHV1-18*02");
}

TEST_F(GermlineDBTest, TestIghChainDbCorrectness) {
    std::string human_ighv_file = "data/germline/human/IG/IGHV.fa";
    std::string human_ighd_file = "data/germline/human/IG/IGHD.fa";
    std::string human_ighj_file = "data/germline/human/IG/IGHJ.fa";
    using namespace germline_utils;
    ChainType chain_type(ImmuneChainType::HeavyIgChain);
    INFO("Creating DB for chain of type " << chain_type);
    ChainDatabase igh_database(chain_type);
    igh_database.AddGenesFromFile(SegmentType::VariableSegment, human_ighv_file);
    igh_database.AddGenesFromFile(SegmentType::DiversitySegment, human_ighd_file);
    igh_database.AddGenesFromFile(SegmentType::JoinSegment, human_ighj_file);
    ASSERT_EQ(igh_database.GenesNumber(SegmentType::VariableSegment), 315);
    INFO("Database contains " << igh_database.GenesNumber(SegmentType::VariableSegment) <<
                 " genes of type " << igh_database.GetDb(SegmentType::VariableSegment).GeneType());
    ASSERT_EQ(igh_database.GenesNumber(SegmentType::DiversitySegment), 34);
    INFO("Database contains " << igh_database.GenesNumber(SegmentType::DiversitySegment) <<
         " genes of type " << igh_database.GetDb(SegmentType::DiversitySegment).GeneType());
    ASSERT_EQ(igh_database.GenesNumber(SegmentType::JoinSegment), 13);
    INFO("Database contains " << igh_database.GenesNumber(SegmentType::JoinSegment) <<
         " genes of type " << igh_database.GetDb(SegmentType::JoinSegment).GeneType());
}

TEST_F(GermlineDBTest, TestTraChainDbCorrectness) {
    std::string human_trav_file = "data/germline/human/TCR/TRAV.fa";
    std::string human_traj_file = "data/germline/human/TCR/TRAJ.fa";
    using namespace germline_utils;
    ChainType chain_type(ImmuneChainType::AlphaTcrChain);
    INFO("Creating DB for chain of type " << chain_type);
    ChainDatabase tra_database(chain_type);
    tra_database.AddGenesFromFile(SegmentType::VariableSegment, human_trav_file);
    tra_database.AddGenesFromFile(SegmentType::JoinSegment, human_traj_file);
    ASSERT_EQ(tra_database.GenesNumber(SegmentType::VariableSegment), 103);
    INFO("Database contains " << tra_database.GenesNumber(SegmentType::VariableSegment) <<
         " genes of type " << tra_database.GetDb(SegmentType::VariableSegment).GeneType());
    ASSERT_EQ(tra_database.GenesNumber(SegmentType::JoinSegment), 68);
    INFO("Database contains " << tra_database.GenesNumber(SegmentType::JoinSegment) <<
         " genes of type " << tra_database.GetDb(SegmentType::JoinSegment).GeneType());
}

TEST_F(GermlineDBTest, TestCustomVariableDbCorrectness) {
    std::string human_ighj_file = "data/germline/human/IG/IGHJ.fa";
    std::string human_igkj_file = "data/germline/human/IG/IGKJ.fa";
    std::string human_iglj_file = "data/germline/human/IG/IGLJ.fa";
    std::string human_traj_file = "data/germline/human/TCR/TRAJ.fa";
    using namespace germline_utils;
    INFO("Creating custom DB for J genes");
    CustomGeneDatabase custom_j_db(SegmentType::JoinSegment);
    custom_j_db.AddDatabase(ImmuneGeneType(ChainType(ImmuneChainType::HeavyIgChain),SegmentType::JoinSegment),
                            human_ighj_file);
    custom_j_db.AddDatabase(ImmuneGeneType(ChainType(ImmuneChainType::LambdaIgChain),SegmentType::JoinSegment),
                            human_igkj_file);
    custom_j_db.AddDatabase(ImmuneGeneType(ChainType(ImmuneChainType::KappaIgChain),SegmentType::JoinSegment),
                            human_iglj_file);
    custom_j_db.AddDatabase(ImmuneGeneType(ChainType(ImmuneChainType::AlphaTcrChain),SegmentType::JoinSegment),
                            human_traj_file);
    ASSERT_EQ(custom_j_db.size(), 100);
    INFO("Custom database of " << custom_j_db.Segment() << " segments contains " << custom_j_db.size() << " records");
}
