#include <gtest/gtest.h>
#include <logger/log_writers.hpp>

#include <sstream>

#include <cdr_config.hpp>
#include <germline_utils/germline_db_generator.hpp>
#include <vj_parallel_processor.hpp>
#include <convert.hpp>

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

vj_finder::VJFinderConfig vj_finder_config;
core::ReadArchive read_archive;
vj_finder::VJAlignmentInfo alignment_info;

class VJFinderTest : public ::testing::Test {
public:
    void SetUp() {
    }
};

void TestVSegmentIgBlastConsistent() {
    INFO("Checking consistency of V hits with IgBlast");
    ASSERT_EQ(core::seqan_string_to_string(alignment_info.GetVJHitsByIndex(0).GetVHitByIndex(0).ImmuneGene().name()),
              "IGHV3-7*01");
    ASSERT_EQ(core::seqan_string_to_string(alignment_info.GetVJHitsByIndex(1).GetVHitByIndex(0).ImmuneGene().name()),
              "IGHV4-34*01");
    ASSERT_EQ(core::seqan_string_to_string(alignment_info.GetVJHitsByIndex(2).GetVHitByIndex(0).ImmuneGene().name()),
              "IGHV1-2*02");
    ASSERT_EQ(core::seqan_string_to_string(alignment_info.GetVJHitsByIndex(3).GetVHitByIndex(0).ImmuneGene().name()),
              "IGHV3-49*04");
    ASSERT_EQ(core::seqan_string_to_string(alignment_info.GetVJHitsByIndex(4).GetVHitByIndex(0).ImmuneGene().name()),
              "IGHV3-30*18");
    ASSERT_EQ(core::seqan_string_to_string(alignment_info.GetVJHitsByIndex(5).GetVHitByIndex(0).ImmuneGene().name()),
              "IGHV1-69*01");
}

void TestJSegmentIgBlastConsistent() {
    INFO("Checking consistency of J hits with IgBlast");
    ASSERT_EQ(core::seqan_string_to_string(alignment_info.GetVJHitsByIndex(0).GetJHitByIndex(0).ImmuneGene().name()),
              "IGHJ4*02");
    ASSERT_EQ(core::seqan_string_to_string(alignment_info.GetVJHitsByIndex(1).GetJHitByIndex(0).ImmuneGene().name()),
              "IGHJ4*02");
    ASSERT_EQ(core::seqan_string_to_string(alignment_info.GetVJHitsByIndex(2).GetJHitByIndex(0).ImmuneGene().name()),
              "IGHJ5*02");
    ASSERT_EQ(core::seqan_string_to_string(alignment_info.GetVJHitsByIndex(3).GetJHitByIndex(0).ImmuneGene().name()),
              "IGHJ6*02");
    ASSERT_EQ(core::seqan_string_to_string(alignment_info.GetVJHitsByIndex(4).GetJHitByIndex(0).ImmuneGene().name()),
              "IGHJ6*02");
    ASSERT_EQ(core::seqan_string_to_string(alignment_info.GetVJHitsByIndex(5).GetJHitByIndex(0).ImmuneGene().name()),
              "IGHJ4*02");
}

void TestStartAndEndAlignmentPositions() {
    INFO("Checking start and end alignment positions");
    for (size_t i = 0; i < alignment_info.NumVJHits(); i++) {
        ASSERT_EQ(alignment_info.GetVJHitsByIndex(i).GetVHitByIndex(0).FirstMatchReadPos(), 0);
        ASSERT_EQ(alignment_info.GetVJHitsByIndex(i).GetJHitByIndex(0).LastMatchReadPos(),
                  read_archive[i].length());
    }
}

void TestReadInversion() {
    INFO("Checking read inversion");
    ASSERT_EQ("GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGTCTGGGGGGTCCCTGAGACTCTCCTGTGCAGC"
              "CTCTGGATTCACCTTTAGTACTTATTGGATGAGCTGGGTCCGCCAGGCTCCAGGGATGGGGCTGGAGTGGG"
              "TGGCCGACATAAAGGAAGATGGAAGTGAGAAATACTATGTGGACTCTGTGAAGGGCCGATTCACCATATCC"
              "AGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTA"
              "CTGTGCGAGAGATTATTGGGGGCCCAATGAGTGGGGCCAGGGAACCACGGTCACCGTCTCCTCAG",
            core::seqan_string_to_string(read_archive[0].seq));
}

void TestReadLeftRightCropping() {
    INFO("Checking left & right cropping");
    ASSERT_EQ("CAGGTGCAGCTACAGCAGTGGGGCGCAGGACTGTTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCGCTGT"
              "CTATGGTGGGTCCTTCAGTGGTTACTACTGGAGCTGGATCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGA"
              "TTGGGGAAATCAATCATAGTGGAAGCACCAACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAGTA"
              "GACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCTGTGTATTACTG"
              "TGCGAGAGGCCTGAGGGCGTTTTCTCGGATACCCCCATTCGAAGACTACTGGGGCCAGGGAACCCTGGTCA"
              "CCGTCTCCTCAG",
              core::seqan_string_to_string(read_archive[1].seq));
}

void TestReadLeftRightFilling() {
    INFO("Checking left & right filling");
    ASSERT_EQ("CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGCGAAGGTCTCCTGCATGGCT"
              "TCTGGATACTCACTCACCGGCTACTATATACACTGGGTGCGACAGGCCCCCGGACAGGGGCTTGAGTGGATGG"
              "GATGGATGAACCCTAACAGCGGTGGCACAACCGATGCACAGAAGTGTCATGGCAGGCTCACCATGACCAGCGA"
                      "CACGTCCATCAACACAGCCTACATGGGGCTGAGCAGGCTGAGATCCGACGACACGGCCGTGTATTA"
                      "CTGTGCGAGAGATCGTACGTATTACGATTTTTGGAGTGCCCCCAGCTCGCGTGGAGGCAACTGGTT"
                      "CGACCCCTGGGGTCAGGGAACCCTGGTCACCGTCTCCTCAG",
              core::seqan_string_to_string(read_archive[2].seq));
}

TEST_F(VJFinderTest, BaseVJFinderTest) {
    create_console_logger();
    std::string config_fname = "configs/vj_finder/config.info";
    vj_finder::load(vj_finder_config, config_fname);
    vj_finder_config.algorithm_params.fix_crop_fill_params.fill_right = true;
    vj_finder_config.algorithm_params.fix_crop_fill_params.fix_right = 3;
    read_archive.ExtractFromFile("test_dataset/vj_finder_test.fastq");
    germline_utils::GermlineDbGenerator db_generator(vj_finder_config.io_params.input_params.germline_input,
                                                     vj_finder_config.algorithm_params.germline_params);
    auto v_gene_database = db_generator.GenerateVariableDb();
    auto j_gene_database = db_generator.GenerateJoinDb();
    vj_finder::VJParallelProcessor processor(read_archive,
                                                vj_finder_config.algorithm_params,
                                                v_gene_database,
                                                j_gene_database,
                                                vj_finder_config.run_params.num_threads);
    alignment_info = processor.Process();

    INFO("Checking number of aligned sequences");
    ASSERT_EQ(alignment_info.NumVJHits(), 6);
    TestVSegmentIgBlastConsistent();
    TestJSegmentIgBlastConsistent();
    TestStartAndEndAlignmentPositions();
    TestReadInversion();
    TestReadLeftRightCropping();
    TestReadLeftRightFilling();
}
