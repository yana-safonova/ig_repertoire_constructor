#include <gtest/gtest.h>
#include <logger/log_writers.hpp>

#include <sstream>

#include <cdr_config.hpp>
#include <germline_db_generator.hpp>
#include <germline_db_labeler.hpp>
#include <vj_parallel_processor.hpp>
#include <read_labeler.hpp>
#include <convert.hpp>

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

cdr_labeler::CDRLabelerConfig config;
annotation_utils::CDRAnnotatedCloneSet annotated_clone_set;

class AntEvoloTest : public ::testing::Test {
public:
    void SetUp() {
        create_console_logger();
        std::string config_fname = "configs/cdr_labeler/config.info";
        config.load(config_fname);
        config.vj_finder_config.algorithm_params.germline_params.loci = "IG";
        vj_finder::GermlineDbGenerator db_generator(config.vj_finder_config.io_params.input_params.germline_input,
                                                    config.vj_finder_config.algorithm_params.germline_params);
        auto v_gene_database = db_generator.GenerateVariableDb();
        auto j_gene_database = db_generator.GenerateJoinDb();
        auto v_labeling = cdr_labeler::GermlineDbLabeler(v_gene_database, config.cdrs_params).ComputeLabeling();
        auto j_labeling = cdr_labeler::GermlineDbLabeler(j_gene_database, config.cdrs_params).ComputeLabeling();
        auto filtered_v_db = v_labeling.CreateFilteredDb();
        auto filtered_j_db = j_labeling.CreateFilteredDb();
        core::ReadArchive read_archive("test_dataset/antevolo_base_test.fastq");
        vj_finder::VJParallelProcessor processor(read_archive, config.vj_finder_config.algorithm_params,
                                                 filtered_v_db, filtered_j_db,
                                                 config.run_params.num_threads);
        auto alignment_info = processor.Process();
        config.shm_params.shm_finding_algorithm =
                cdr_labeler::CDRLabelerConfig::SHMFindingParams::SHMFindingAlgorithm::CDRFilteringSHMAlgorithm;
        cdr_labeler::ReadCDRLabeler read_labeler(config.shm_params, v_labeling, j_labeling);
        annotated_clone_set = read_labeler.CreateAnnotatedCloneSet(alignment_info);
    }
};

TEST_F(AntEvoloTest, AntEvoloBaseTest) {
    ASSERT_EQ(annotated_clone_set[0].VSHMs().size(), 11);
}