//
// Created by Andrew Bzikadze on 7/11/16.
//

#include <gtest/gtest.h>
#include "vdj_config.hpp"
#include <read_archive.hpp>
#include <germline_db_generator.hpp>
#include <vj_parallel_processor.hpp>
#include <block_alignment/block_alignment_converter.hpp>
#include <vdj_alignments/vdj_hits.hpp>
#include <vdj_alignments/vdj_hits_storage.hpp>
#include <alignment_utils/alignment_positions.hpp>
#include <vdj_alignments/aligners/simple_d_aligner.hpp>
#include <model/recombination_model.hpp>
#include <germline_utils/chain_type.hpp>
#include <recombination_calculator/hc_model_based_recombination_calculator.hpp>
#include "vdj_alignments/aligners/simple_d_aligner.hpp"

using namespace vdj_labeler;
// using namespace recombination_utils;
// using namespace alignment_utils;
using namespace germline_utils;
using namespace seqan;
using namespace core;
using namespace std;


class AlignersTest: public ::testing::Test {
public:
    core::ReadArchive read_archive;
    CustomGeneDatabase v_db;
    CustomGeneDatabase d_db;
    CustomGeneDatabase j_db;

    AlignersTest() :
        v_db(SegmentType::VariableSegment),
        d_db(SegmentType::DiversitySegment),
        j_db(SegmentType::JoinSegment)
    { }

    void SetUp() {
        std::string cfg_filename = "configs/vdj_labeler_0.3/configs.info";
        vdj_labeler::VDJLabelerConfig config;
        config.load(cfg_filename);
        std::string input_filename = config.io_params.input_params.input_sequences;
        std::string v_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.variable_genes;
        std::string d_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.diversity_genes;
        std::string j_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.join_genes;

        read_archive = core::ReadArchive(input_filename);
        read_archive.FixSpacesInHeaders();

        v_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::VariableSegment), v_germline_genes_fname);
        d_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::DiversitySegment), d_germline_genes_fname);
        j_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::JoinSegment), j_germline_genes_fname);
    }
};

TEST_F(AlignersTest, SimpleDAlignerTest) {
    for (auto read_it = read_archive.cbegin(); read_it != read_archive.cend(); ++read_it) {
        for (size_t i = 0; i < d_db.size(); ++i) {
            alignment_utils::AlignmentPositions alignment_positions(
                std::make_pair<size_t, size_t>(0, read_it->length() - 1),
                std::make_pair<size_t, size_t>(0, 10)); // 0 and 10 don't play any role

            alignment_utils::ImmuneGeneAlignmentPositions immune_alignment_positions(alignment_positions,
                                                                                     d_db[i],
                                                                                     *read_it);
            auto alignment = SimpleDAligner().ComputeAlignment(immune_alignment_positions);
            EXPECT_LE(immune_alignment_positions.ReadStartPos(), alignment.StartQueryPosition());
            EXPECT_GE(immune_alignment_positions.ReadEndPos(), alignment.EndQueryPosition() - alignment.NumberGaps());
            EXPECT_EQ(seqan::length(seqan::row(alignment.Alignment(), 0)),
                      seqan::length(seqan::row(alignment.Alignment(), 1)));
            EXPECT_EQ(seqan::length(seqan::row(alignment.Alignment(), 0)), alignment.SubjectAlignmentLength());
            EXPECT_EQ(seqan::length(seqan::row(alignment.Alignment(), 1)), alignment.QueryAlignmentLength());
        }
    }
}
