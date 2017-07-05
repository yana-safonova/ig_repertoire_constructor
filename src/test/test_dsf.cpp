//
// Created by Yana Safonova on 14/02/16.
//

#include <gtest/gtest.h>
//#include <boost/filesystem.hpp>

#include <logger/log_writers.hpp>
#include "../dense_sgraph_finder/graph_decomposer/dense_subgraph_constructor.hpp"
#include "../dense_sgraph_finder/dsf_config.hpp"
#include "../graph_utils/graph_io.hpp"

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

dsf_config::dense_sgraph_finder_params CreateStandardDsfParams() {
    dsf_config::dense_sgraph_finder_params dsf_params;
    dsf_params.min_supernode_size = 5;
    dsf_params.create_trivial_decomposition = false;
    dsf_params.min_fillin_threshold = 0.6;
    dsf_params.min_graph_size = 5;
    dsf_params.primary_edge_fillin = 0.3;
    return dsf_params;
}

dsf_config::metis_io_params CreateStandardMetisParams(string output_dir) {
    dsf_config::metis_io_params metis_params;
    metis_params.path_to_metis = "build/release/bin/";
    metis_params.run_metis = path::append_path(metis_params.path_to_metis, "./metis");
    metis_params.trash_output = path::append_path(output_dir, "metis.output");
    return metis_params;
}

SparseGraphPtr test_graph;
std::string output_dir;

class DsfTest : public ::testing::Test {
public:
    void SetUp();

    void TearDown();
};

void DsfTest::SetUp() {
    create_console_logger();
    output_dir = "dsf_unit_tests";
    // Yana: boost_filesystem does not compile on mac
    path::make_dir(output_dir);
    // todo: add graph with weighted vertices in test
    std::string graph_filename = "test_dataset/dsf/weighted_vertices.graph";
    GraphReader graph_reader(graph_filename);
    test_graph = graph_reader.CreateGraph();
}

void DsfTest::TearDown() {
    // Yana: boost_filesystem does not compile on mac
    path::remove_dir(output_dir);
}

// create decomposition into dense subgraphs
// check that each constructed dense subgraph contains at most one supernode
TEST_F(DsfTest, TestSupernodesAreSeparated) {
    dsf_config::dense_sgraph_finder_params dsf_params = CreateStandardDsfParams();
    auto metis_io_params = CreateStandardMetisParams(output_dir);
    dense_subgraph_finder::MetisDenseSubgraphConstructor dsf_constructor(dsf_params,
                                                                         metis_io_params,
                                                                         path::append_path(output_dir,
                                                                                           "graph_copy.graph"));
    DecompositionPtr decomposition = dsf_constructor.CreateDecomposition(test_graph);
    INFO(decomposition->Size() << " dense subgraphs were constructed");
    for(size_t i = 0; i < decomposition->Size(); i++) {
        auto cur_class = decomposition->GetClass(i);
        size_t num_supernodes = 0;
        for(auto v = cur_class.begin(); v != cur_class.end(); v++)
            if(test_graph->WeightOfVertex(*v) >= dsf_params.min_supernode_size)
                num_supernodes += 1;
        ASSERT_EQ(num_supernodes <= 1, true);
    }
    INFO("Each dense subgraph contains at most one supernode");
}
