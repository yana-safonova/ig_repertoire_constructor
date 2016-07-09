//
// Created by Andrew Bzikadze on 7/9/16.
//

#include <gtest/gtest.h>
#include "vdj_config.hpp"
#include "model/recombination_model.hpp"
#include "read_archive.hpp"
#include "vdj_alignments/vdj_hits.hpp"
#include "recombination_utils/cleaved_gene.hpp"

#include <string>

using namespace vdj_labeler;
// using namespace recombination_utils;
// using namespace alignment_utils;
using namespace germline_utils;
using namespace seqan;
using namespace core;
using namespace std;

germline_utils::ChainDatabase hc_db(germline_utils::ImmuneChainType::HeavyIgChain);
std::string input_filename;

class RecombinationModelTest : public ::testing::Test {
    void SetUp() {
        std::string cfg_filename = "configs/vdj_labeler_0.3/configs.info";
        vdj_labeler::VDJLabelerConfig config;
        config.load(cfg_filename);
        std::string input_filename = config.io_params.input_params.input_sequences;
        std::string v_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.variable_genes;
        std::string d_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.diversity_genes;
        std::string j_germline_genes_fname = config.io_params.input_params.germlines.igh_genes.join_genes;

        hc_db.AddGenesFromFile(SegmentType::VariableSegment, v_germline_genes_fname);
        hc_db.AddGenesFromFile(SegmentType::DiversitySegment, d_germline_genes_fname);
        hc_db.AddGenesFromFile(SegmentType::JoinSegment, j_germline_genes_fname);
    }
};

TEST_F(RecombinationModelTest, TestIgGeneProbabilityModel) {
    std::ifstream in("src/vdj_labeler_0.3/test/blank_model.csv");
    HCProbabilityRecombinationModel model(in, hc_db);
    ASSERT_EQ(model.GetVGeneProbabilityModel().size(), 97);
    EXPECT_DOUBLE_EQ(model.GetVGeneProbabilityModel().GetProbabilityByGenId(0), 0.036479279395);
    EXPECT_DOUBLE_EQ(model.GetVGeneProbabilityModel().GetProbabilityByGenId(9), 0.021516909466);

    ASSERT_EQ(model.GetDGeneProbabilityModel().size(), 35);
    EXPECT_DOUBLE_EQ(model.GetDGeneProbabilityModel().GetProbabilityByGenId(4), 0.004103032);
    EXPECT_DOUBLE_EQ(model.GetProbabilityByGenId(SegmentType::DiversitySegment, 4), 0.004103032);

    ASSERT_EQ(model.GetJGeneProbabilityModel().size(), 7);
    EXPECT_DOUBLE_EQ(model.GetProbabilityByGenId(SegmentType::JoinSegment, 4), 0.13860616);

    ASSERT_EQ(model.GetVDNongenomicInsertionModel().GetInsertionProbabilities().size(), 61);
    EXPECT_DOUBLE_EQ(model.GetVDNongenomicInsertionModel().GetInsertionProbabilityByLength(60), 0.000100717807129166);
    EXPECT_DOUBLE_EQ(model.GetVDNongenomicInsertionModel().GetInsertionProbabilityByLength(29), 0.00402026064657394);

    ASSERT_EQ(model.GetVDNongenomicInsertionModel().GetTransitionMatrix().size(), 4);
    EXPECT_DOUBLE_EQ(model.GetVDNongenomicInsertionModel().GetTransitionProbability('A', 'C'), 0.158188525976022);
    EXPECT_DOUBLE_EQ(model.GetVDNongenomicInsertionModel().GetTransitionProbability('A', 'G'), 0.220189234667469);
    EXPECT_DOUBLE_EQ(model.GetVDNongenomicInsertionModel().GetTransitionProbability(make_pair('A', 'G')), 0.220189234667469);

    ASSERT_EQ(model.GetDJNongenomicInsertionModel().GetInsertionProbabilities().size(), 61);
    EXPECT_DOUBLE_EQ(model.GetDJNongenomicInsertionModel().GetInsertionProbabilityByLength(60), 0.000212678981597781);
    EXPECT_DOUBLE_EQ(model.GetDJNongenomicInsertionModel().GetInsertionProbabilityByLength(44), 0.00152129688100844);

    ASSERT_EQ(model.GetDJNongenomicInsertionModel().GetTransitionMatrix().size(), 4);
    EXPECT_DOUBLE_EQ(model.GetDJNongenomicInsertionModel().GetTransitionProbability('C', 'C'), 0.388605972280164);

    ASSERT_EQ(model.GetVPalindromeDeletionModel().size(), 97);
    EXPECT_EQ(model.GetVPalindromeDeletionModel().GetDeletionLength().front(), -6);
    EXPECT_EQ(model.GetVPalindromeDeletionModel().GetDeletionLength().back(), 40);
    EXPECT_EQ(model.GetVPalindromeDeletionModel().GetDeletionProbability(0, -6), 0.0000467021538387012);
    EXPECT_EQ(model.GetVPalindromeDeletionModel().GetDeletionProbability(2, -4), 0.00271564283179921);

    ASSERT_EQ(model.GetDLeftPalindromeDeletionModel().size(), 35);
    EXPECT_EQ(model.GetDLeftPalindromeDeletionModel().GetDeletionLength().front(), -6);
    EXPECT_EQ(model.GetDLeftPalindromeDeletionModel().GetDeletionLength().back(), 37);
    EXPECT_EQ(model.GetDLeftPalindromeDeletionModel().GetDeletionProbability(0, -6), 0.000435083);
    EXPECT_EQ(model.GetDLeftPalindromeDeletionModel().GetDeletionProbability(2, -4), 0.00083265410);

    ASSERT_EQ(model.GetDRightPalindromeDeletionModel().size(), 35);
    EXPECT_EQ(model.GetDRightPalindromeDeletionModel().GetDeletionLength().front(), -6);
    EXPECT_EQ(model.GetDRightPalindromeDeletionModel().GetDeletionLength().back(), 37);
    EXPECT_EQ(model.GetDRightPalindromeDeletionModel().GetDeletionProbability(0, -6), 0.000118509);
    EXPECT_EQ(model.GetDRightPalindromeDeletionModel().GetDeletionProbability(2, -4), 0.00041302831);

    ASSERT_EQ(model.GetJPalindromeDeletionModel().size(), 7);
    EXPECT_EQ(model.GetJPalindromeDeletionModel().GetDeletionLength().front(), -6);
    EXPECT_EQ(model.GetJPalindromeDeletionModel().GetDeletionLength().back(), 41);
    EXPECT_EQ(model.GetJPalindromeDeletionModel().GetDeletionProbability(0, -6), 0);
    EXPECT_EQ(model.GetJPalindromeDeletionModel().GetDeletionProbability(2, -4), 0.00320774662353228);

    // HCModelBasedRecombinationCalculator recombination_calculator(model);
}

TEST_F(RecombinationModelTest, TestIgGeneProbabilityModelCalculator) {
    // TODO complete this test when the rest of the codebase is ready.
    // std::ifstream in("src/vdj_labeler_0.3/test/blank_model.csv");
    // HCProbabilityRecombinationModel model(in, hc_db);

    // core::ReadArchive reads_archive(input_filename);
    // size_t read_index = 3;
    // ReadPtr read = reads_archive[read_index];
    // VDJHitsPtrhits_3 = (*hits_storage)[read_index];
    // INFO("Read 3. #V: " << hits_3->VHitsNumber() <<
    //     ", #D: " << hits_3->DHitsNumber() <<
    //     ", #J: " << hits_3->JHitsNumber());

    // auto v_alignment = hits_3->GetAlignmentByIndex(IgGeneType::variable_gene, 0);
    // CleavedIgGeneAlignment v_event_0(v_alignment, 0, 0, 0, 0);
    // CleavedIgGeneAlignment v_event_1(v_alignment, 0, -1, 0, 0);
    // CleavedIgGeneAlignment v_event_2(v_alignment, 0, -2, 0, 0);
    // CleavedIgGeneAlignment v_event_3(v_alignment, 0, -3, 0, 1);

    // auto d_alignment = hits_3->GetAlignmentByIndex(IgGeneType::diversity_gene, 0);
    // CleavedIgGeneAlignment d_event_0(d_alignment, 1, 8, 0, 0);

    // auto j_alignment = hits_3->GetAlignmentByIndex(IgGeneType::join_gene, 0);
    // CleavedIgGeneAlignment j_event_0(j_alignment, 0, 0, 1, 0);
    // CleavedIgGeneAlignment j_event_1(j_alignment, 1, 0, 0, 0);

    // NongenomicInsertion vd_insertion_0(425, 441);
    // NongenomicInsertion vd_insertion_1(426, 441);
    // NongenomicInsertion vd_insertion_2(427, 441);
    // NongenomicInsertion vd_insertion_3(428, 441);

    // NongenomicInsertion dj_insertion_0(453, 452);
    // NongenomicInsertion dj_insertion_1(453, 453);
    // RecombinationStorage<HCRecombination> recombination_storage(read_3);
    // recombination_storage.AddRecombination(HCRecombination(read_3, v_event_0, d_event_0, j_event_0,
    //                                                        vd_insertion_0, dj_insertion_0));
    // recombination_storage.AddRecombination(HCRecombination(read_3, v_event_1, d_event_0, j_event_0,
    //                                                        vd_insertion_1, dj_insertion_0));
    // recombination_storage.AddRecombination(HCRecombination(read_3, v_event_2, d_event_0, j_event_0,
    //                                                        vd_insertion_2, dj_insertion_0));
    // recombination_storage.AddRecombination(HCRecombination(read_3, v_event_3, d_event_0, j_event_0,
    //                                                        vd_insertion_3, dj_insertion_0));
    // recombination_storage.AddRecombination(HCRecombination(read_3, v_event_0, d_event_0, j_event_1,
    //                                                        vd_insertion_0, dj_insertion_1));
    // recombination_storage.AddRecombination(HCRecombination(read_3, v_event_1, d_event_0, j_event_1,
    //                                                        vd_insertion_1, dj_insertion_1));
    // recombination_storage.AddRecombination(HCRecombination(read_3, v_event_2, d_event_0, j_event_1,
    //                                                        vd_insertion_2, dj_insertion_1));
    // recombination_storage.AddRecombination(HCRecombination(read_3, v_event_3, d_event_0, j_event_1,
    //                                                        vd_insertion_3, dj_insertion_1));
    // INFO(recombination_storage.size() << " recombinaions were generated");
}
