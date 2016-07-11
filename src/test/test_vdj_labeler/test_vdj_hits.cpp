//
// Created by Andrew Bzikadze on 7/9/16.
//

#include <gtest/gtest.h>
#include "vdj_config.hpp"
#include "model/recombination_model.hpp"
#include "read_archive.hpp"
#include "vdj_alignments/vdj_hits.hpp"
#include "recombination_utils/cleaved_gene.hpp"
#include <vj_parallel_processor.hpp>
#include "vdj_alignments/vdj_hits_storage.hpp"

using namespace vdj_labeler;
// using namespace recombination_utils;
// using namespace alignment_utils;
using namespace germline_utils;
using namespace seqan;
using namespace core;
using namespace std;

vj_finder::VJAlignmentInfo alignment_info;

class VDJHitsTest: public ::testing::Test {
public:
    void SetUp() {
        std::string cfg_filename = "configs/vdj_labeler_0.3/configs.info";
        vdj_labeler::VDJLabelerConfig config;
        config.load(cfg_filename);
        std::string input_filename = config.io_params.input_params.input_sequences;
        std::string v_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.variable_genes;
        std::string d_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.diversity_genes;
        std::string j_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.join_genes;

        core::ReadArchive read_archive(input_filename);
        read_archive.FixSpacesInHeaders();

        CustomGeneDatabase v_db(SegmentType::VariableSegment);
        CustomGeneDatabase d_db(SegmentType::DiversitySegment);
        CustomGeneDatabase j_db(SegmentType::JoinSegment);

        v_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::VariableSegment), v_germline_genes_fname);
        d_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::DiversitySegment), d_germline_genes_fname);
        j_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::JoinSegment), j_germline_genes_fname);

        germline_utils::ChainDatabase hc_db(germline_utils::ImmuneChainType::HeavyIgChain);
        hc_db.AddGenesFromFile(SegmentType::VariableSegment, v_germline_genes_fname);
        hc_db.AddGenesFromFile(SegmentType::DiversitySegment, d_germline_genes_fname);
        hc_db.AddGenesFromFile(SegmentType::JoinSegment, j_germline_genes_fname);

        vj_finder::VJParallelProcessor processor(read_archive, config.vj_finder_config.algorithm_params,
                                                 v_db, j_db,
                                                 config.run_params.threads_count);
        alignment_info = processor.Process();
    }
};


TEST_F(VDJHitsTest, VDJHitsBaseTest) {
    auto vdj_hits = VDJHits(alignment_info.AlignmentRecords()[0]);
    EXPECT_EQ(vdj_hits.VHits().GeneType(), SegmentType::VariableSegment);
    EXPECT_EQ(vdj_hits.JHits().GeneType(), SegmentType::JoinSegment);

    EXPECT_EQ(vdj_hits.Read()->name, "MIG_UMI:ACAGGCTAGAGAAC:10_25811");
}

TEST_F(VDJHitsTest, VJDJHitsStorageTest) {
    auto vdj_storage = VDJHitsStorage(alignment_info);
    EXPECT_EQ(vdj_storage.size(), alignment_info.AlignmentRecords().size());
    for(auto it = vdj_storage.cbegin(); it != vdj_storage.cend(); ++it) {
        EXPECT_EQ((*it)->VHits().GeneType(), SegmentType::VariableSegment);
        EXPECT_EQ((*it)->JHits().GeneType(), SegmentType::JoinSegment);
    }
}
