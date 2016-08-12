//
// Created by Andrew Bzikadze on 8/10/16.
//

#include <logger/log_writers.hpp>
#include <gtest/gtest.h>
#include "../../vdj_labeler_0.3/vdj_config.hpp"
#include "vdj_config.hpp"
#include <read_archive.hpp>
#include <germline_db_generator.hpp>
#include <vj_parallel_processor.hpp>
#include <block_alignment/block_alignment_converter.hpp>
#include <vdj_alignments/vdj_hits/vdj_hits.hpp>
#include <vdj_alignments/vdj_hits/vdj_hits_storage.hpp>
#include <alignment_utils/alignment_positions.hpp>
#include <germline_utils/chain_type.hpp>
#include "vdj_alignments/aligners/simple_d_aligner.hpp"
#include "vdj_alignments/hits_calculator/alignment_quality_checkers/match_threshold_alignment_quality_checker.hpp"
#include "vdj_alignments/hits_calculator/d_alignment_positions_checkers/info_based_d_alignment_position_checker.hpp"
#include "vdj_alignments/hits_calculator/d_alignment_positions_calculator/custom_d_alignment_positions_calculator.hpp"
#include "vdj_alignments/hits_calculator/d_hits_calculator/info_based_d_hits_calculator.hpp"

using namespace vdj_labeler;
// using namespace recombination_utils;
// using namespace alignment_utils;
using namespace germline_utils;
using namespace seqan;
using namespace core;
using namespace std;

class DHitsCalculatorTest : public ::testing::Test {
    void SetUp() {
        using namespace logging;
        logger *lg = create_logger("");
        lg->add_writer(std::make_shared<console_writer>());
        attach_logger(lg);
    }
};

TEST_F(DHitsCalculatorTest, DynamicsTest) {
    std::string cfg_filename = "configs/vdj_labeler_0.3/configs.info";
    vdj_labeler::VDJLabelerConfig config;
    config.load(cfg_filename);
    std::string input_filename = "test_dataset/vdj_labeler_0.3/LongCDR3.fasta";
    std::string v_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.variable_genes;
    std::string d_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.diversity_genes;
    std::string j_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.join_genes;

    auto read_archive = core::ReadArchive(input_filename);
    read_archive.FixSpacesInHeaders();

    CustomGeneDatabase v_db(SegmentType::VariableSegment);
    CustomGeneDatabase d_db(SegmentType::DiversitySegment);
    CustomGeneDatabase j_db(SegmentType::JoinSegment);

    v_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::VariableSegment), v_germline_genes_fname);
    d_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::DiversitySegment), d_germline_genes_fname);
    j_db.AddDatabase(ImmuneGeneType(ChainType("IGH"), SegmentType::JoinSegment), j_germline_genes_fname);

    vj_finder::VJParallelProcessor processor(read_archive, config.vj_finder_config.algorithm_params,
                                             v_db, j_db,
                                             config.run_params.threads_count);
    vj_finder::VJAlignmentInfo alignment_info = processor.Process();

    SimpleDAligner d_aligner;
    MatchThresholdAlignmentQualityChecker quality_checker(5);
    // ThresholdAlignmentQualityChecker quality_checker(1);
    InfoBasedDAlignmentPositionChecker position_checker(config.d_align_quality_params);
    CustomDAlignmentPositionsCalculator positions_calculator;
    InfoBasedDHitsCalculator calculator(
        d_db.GetConstDbByGeneType(ImmuneGeneType(ChainType("IGH"), SegmentType::DiversitySegment)),
        d_aligner, quality_checker, position_checker, positions_calculator, config.d_align_quality_params);

    VDJHitsStorage vdj_storage2(alignment_info, calculator);
    EXPECT_EQ(vdj_storage2.size(), 1);
    auto& vdj_hit = vdj_storage2[0];
    EXPECT_EQ(vdj_hit.DHitsNumber(), 1);
    auto& d_hit = vdj_hit.DHits()[0];
    auto& vector_d_alignments = d_hit.DGenesHits();
    EXPECT_TRUE(isEqual(vector_d_alignments[0].Subject().name(), "IGHD6-19*01"));
    EXPECT_TRUE(isEqual(vector_d_alignments[1].Subject().name(), "IGHD3-3*01"));
    EXPECT_TRUE(isEqual(vector_d_alignments[2].Subject().name(), "IGHD1-26*01"));

    EXPECT_EQ(vector_d_alignments[0].AlignmentLength(), 10);
    EXPECT_EQ(vector_d_alignments[1].AlignmentLength(), 12);
    EXPECT_EQ(vector_d_alignments[2].AlignmentLength(), 7);
}