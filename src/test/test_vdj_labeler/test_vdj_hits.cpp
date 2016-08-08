//
// Created by Andrew Bzikadze on 7/9/16.
//

#include <gtest/gtest.h>
#include "vdj_config.hpp"
#include "model/recombination_model.hpp"
#include "read_archive.hpp"
#include "vdj_alignments/vdj_hits/vdj_hits.hpp"
#include "recombination_utils/cleaved_gene.hpp"
#include "vj_parallel_processor.hpp"
#include "vdj_alignments/vdj_hits/vdj_hits_storage.hpp"
#include "vdj_alignments/aligners/simple_d_aligner.hpp"
#include "vdj_alignments/hits_calculator/d_hits_calculator/info_based_d_hits_calculator.hpp"
#include "vdj_alignments/hits_calculator/alignment_quality_checkers/match_threshold_alignment_quality_checker.hpp"
#include "vdj_alignments/hits_calculator/d_alignment_positions_calculator/custom_d_alignment_positions_calculator.hpp"
#include "vdj_alignments/hits_calculator/d_alignment_positions_checkers/info_based_d_alignment_position_checker.hpp"

using namespace vdj_labeler;
// using namespace recombination_utils;
// using namespace alignment_utils;
using namespace germline_utils;
using namespace seqan;
using namespace core;
using namespace std;

vj_finder::VJAlignmentInfo alignment_info;
vdj_labeler::VDJLabelerConfig config;
CustomGeneDatabase v_db(SegmentType::VariableSegment);
CustomGeneDatabase d_db(SegmentType::DiversitySegment);
CustomGeneDatabase j_db(SegmentType::JoinSegment);

class VDJHitsTest: public ::testing::Test {
public:
    void SetUp() {
        std::string cfg_filename = "configs/vdj_labeler_0.3/configs.info";
        config.load(cfg_filename);
        std::string input_filename = config.io_params.input_params.input_sequences;
        std::string v_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.variable_genes;
        std::string d_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.diversity_genes;
        std::string j_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.join_genes;

        core::ReadArchive read_archive(input_filename);
        read_archive.FixSpacesInHeaders();

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
}

TEST_F(VDJHitsTest, VJDJHitsStorageTest) {
    SimpleDAligner d_aligner;
    MatchThresholdAlignmentQualityChecker quality_checker;
    InfoBasedDAlignmentPositionChecker position_checker(config.d_align_quality_params);
    CustomDAlignmentPositionsCalculator positions_calculator;
    InfoBasedDHitsCalculator calculator(
        d_db.GetConstDbByGeneType(ImmuneGeneType(ChainType("IGH"), SegmentType::DiversitySegment)),
        d_aligner, quality_checker, position_checker, positions_calculator);

    VDJHitsStorage vdj_storage (alignment_info);
    // VDJHitsStorage vdj_storage2(alignment_info, calculator);

    EXPECT_EQ(vdj_storage.size(), alignment_info.AlignmentRecords().size());
    for (const auto& vdj_hit : vdj_storage) {
        EXPECT_EQ(vdj_hit.VHits().GeneType(), SegmentType::VariableSegment);
        // EXPECT_EQ(vdj_hit.DHits().GeneType(), SegmentType::DiversitySegment);
        EXPECT_EQ(vdj_hit.JHits().GeneType(), SegmentType::JoinSegment);
    }

    // ASSERT_EQ(vdj_storage.size(), vdj_storage2.size());
    for (size_t hit_n = 0; hit_n < vdj_storage.size(); ++hit_n) {
        auto &vdj_hits = vdj_storage[hit_n];
        // auto &vdj_hits2 = vdj_storage2[hit_n];
        for (const auto &v_hit : vdj_hits.VHits()) {
            EXPECT_EQ(v_hit.Subject().Segment(), SegmentType::VariableSegment);
        }
        // for (const auto &d_hit : vdj_hits.DHits()) {
        //     EXPECT_EQ(d_hit.Subject().Segment(), SegmentType::DiversitySegment);
        // }
        for (const auto &j_hit : vdj_hits.JHits()) {
            EXPECT_EQ(j_hit.Subject().Segment(), SegmentType::JoinSegment);
        }
        // for (const auto &v_hit : vdj_hits2.VHits()) {
        //     EXPECT_EQ(v_hit.Subject().Segment(), SegmentType::VariableSegment);
        // }
        // for (const auto &d_hit : vdj_hits2.DHits()) {
        //     EXPECT_EQ(d_hit.Subject().Segment(), SegmentType::DiversitySegment);
        // }
        // for (const auto &j_hit : vdj_hits2.JHits()) {
        //     EXPECT_EQ(j_hit.Subject().Segment(), SegmentType::JoinSegment);
        // }

        for (size_t i = 0; i < vdj_hits.VHitsNumber(); ++i) {
            EXPECT_TRUE(isEqual(vdj_hits.VHits()[i].Subject().seq(),
                                alignment_info.AlignmentRecords()[hit_n].VHits()[i].ImmuneGene().seq()));
            // EXPECT_TRUE(isEqual(vdj_hits2.VHits()[i].Subject().seq(),
            //                     alignment_info.AlignmentRecords()[hit_n].VHits()[i].ImmuneGene().seq()));
        }
        for (size_t i = 0; i < vdj_hits.JHitsNumber(); ++i) {
            EXPECT_TRUE(isEqual(vdj_hits.JHits()[i].Subject().seq(),
                                alignment_info.AlignmentRecords()[hit_n].JHits()[i].ImmuneGene().seq()));
            // EXPECT_TRUE(isEqual(vdj_hits2.JHits()[i].Subject().seq(),
            //                     alignment_info.AlignmentRecords()[hit_n].JHits()[i].ImmuneGene().seq()));
        }
        // ASSERT_EQ(vdj_hits.VHitsNumber(), vdj_hits2.VHitsNumber());
        for (size_t i = 0; i < vdj_hits.VHitsNumber(); ++i) {
            EXPECT_TRUE(isEqual(vdj_hits.VHits()[i].Query().seq,
                                alignment_info.AlignmentRecords()[hit_n].VHits()[i].Read().seq));
            // EXPECT_TRUE(isEqual(vdj_hits2.VHits()[i].Query().seq,
            //                     alignment_info.AlignmentRecords()[hit_n].VHits()[i].Read().seq));
        }
        // ASSERT_EQ(vdj_hits.JHitsNumber(), vdj_hits2.JHitsNumber());
        for (size_t i = 0; i < vdj_hits.JHitsNumber(); ++i) {
            EXPECT_TRUE(isEqual(vdj_hits.JHits()[i].Query().seq,
                                alignment_info.AlignmentRecords()[hit_n].JHits()[i].Read().seq));
            // EXPECT_TRUE(isEqual(vdj_hits2.JHits()[i].Query().seq,
            //                     alignment_info.AlignmentRecords()[hit_n].JHits()[i].Read().seq));
        }

        // for (size_t i = 0; i < vdj_hits2.DHitsNumber(); ++i) {
        //     EXPECT_TRUE(isEqual(vdj_hits2.DHits()[i].Query().seq,
        //                         alignment_info.AlignmentRecords()[hit_n].VHits()[0].Read().seq));
        // }
    }
}
