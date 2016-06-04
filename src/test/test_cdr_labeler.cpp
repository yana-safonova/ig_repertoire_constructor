#include <gtest/gtest.h>
#include <logger/log_writers.hpp>

#include <sstream>

#include <cdr_config.hpp>
#include <germline_db_generator.hpp>
#include <germline_db_labeler.hpp>

void create_console_logger() {
        using namespace logging;
        logger *lg = create_logger("");
        lg->add_writer(std::make_shared<console_writer>());
        attach_logger(lg);
}

cdr_labeler::CDRLabelerConfig config;
germline_utils::CustomGeneDatabase v_gene_database(germline_utils::SegmentType::VariableSegment);

class CDRLabelerTest : public ::testing::Test {
public:
    void SetUp() {
        create_console_logger();
        std::string config_fname = "configs/cdr_labeler/config.info";
        config.load(config_fname);
        vj_finder::GermlineDbGenerator db_generator(config.vj_finder_config.io_params.input_params.germline_input,
                                                    config.vj_finder_config.algorithm_params.germline_params);
        v_gene_database = db_generator.GenerateVariableDb();
    }
};

std::string GetCDRSequence(const germline_utils::ImmuneGene &immune_gene, annotation_utils::CDRRange cdr_range) {
    auto seqan_cdr = seqan::infixWithLength(immune_gene.seq(), cdr_range.start_pos, cdr_range.length());
    std::stringstream ss;
    ss << seqan_cdr;
    return ss.str();
}

TEST_F(CDRLabelerTest, GermlineCDRsAreConsistentWithIgBlast) {
    auto v_labeling = cdr_labeler::GermlineDbLabeler(v_gene_database, config.cdrs_params).ComputeLabeling();
    auto filtered_v_db = v_labeling.CreateFilteredDb();
    auto ighv1_18 = filtered_v_db[0];
    INFO("Testing CDR1 and CDR2 for gene " << ighv1_18.name());
    auto v_cdrs = v_labeling.GetLabelingByGene(ighv1_18);
    std::string computed_cdr1 = GetCDRSequence(ighv1_18, v_cdrs.cdr1);
    std::string computed_cdr2 = GetCDRSequence(ighv1_18, v_cdrs.cdr2);
    std::string igblast_cdr1 = "GGTTACACCTTTACCAGCTATGGT";
    std::string igblast_cdr2 = "ATCAGCGCTTACAATGGTAACACA";
    ASSERT_EQ(igblast_cdr1, computed_cdr1);
    ASSERT_EQ(igblast_cdr2, computed_cdr2);
}
