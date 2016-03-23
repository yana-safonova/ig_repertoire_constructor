//
// Created by Yana Safonova on 23/03/16.
//

#include <gtest/gtest.h>
#include <logger/log_writers.hpp>

#include "../vdj_utils/germline_utils/germline_databases/immune_gene_database.hpp"
#include "../vdj_utils/germline_utils/germline_databases/chain_database.hpp"
#include "../vdj_utils/germline_utils/germline_databases/custom_gene_database.hpp"

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

class GermlineDBTest : public ::testing::Test {
public:
    void SetUp();

    void TearDown() { }
};

void GermlineDBTest::SetUp() {
    create_console_logger();
    //output_dir = "germline_db_unit_tests";
    //path::make_dir(output_dir);
    // todo: add graph with weighted vertices in test
    //std::string graph_filename = "test_dataset/dsf/weighted_vertices.graph";
    //GraphReader graph_reader(graph_filename);
    //test_graph = graph_reader.CreateGraph();
}

// create decomposition into dense subgraphs
// check that each constructed dense subgraph contains at most one supernode
TEST_F(GermlineDBTest, TestImmuneGeneDbCorrectness) {
    std::string human_ighv_file = "data/germline/human/IG/IGHV.fa";
    using namespace germline_utils;
    ImmuneGeneDatabase ighv_database(ImmuneGeneType(
            ChainType(ImmuneChainType::HeavyIgChain),
            SegmentType::VariableSegment));
    ighv_database.AddGenesFromFile(human_ighv_file);
    INFO("Reading database from " << human_ighv_file);
    ASSERT_EQ(ighv_database.size(), 351);
    INFO("DB contains " << ighv_database.size() << " records");
    std::stringstream chain_str;
    chain_str << ighv_database.Chain();
    ASSERT_EQ(chain_str.str(), "IGH");
    INFO("DB chain type: " << ighv_database.Chain());
    std::stringstream gene_str;
    gene_str << ighv_database.GeneType();
    ASSERT_EQ(gene_str.str(), "IGHV");
    INFO("DB gene type: " << ighv_database.GeneType());
}

TestF(GermlineDBTest, TestChainDbCorrectness) {
    
}

