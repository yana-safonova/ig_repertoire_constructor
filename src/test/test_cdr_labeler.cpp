#include <gtest/gtest.h>
#include <logger/log_writers.hpp>

#include <sstream>

#include <cdr_config.hpp>
#include <germline_utils/germline_db_generator.hpp>
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
germline_utils::CustomGeneDatabase filtered_v_db(germline_utils::SegmentType::VariableSegment);
germline_utils::CustomGeneDatabase filtered_j_db(germline_utils::SegmentType::JoinSegment);

class CDRLabelerTest : public ::testing::Test {
public:
    void SetUp() {
        create_console_logger();
        std::string config_fname = "configs/cdr_labeler/config.info";
        config.load(config_fname);
        config.vj_finder_config.algorithm_params.germline_params.loci = "IG";
        germline_utils::GermlineDbGenerator db_generator(config.vj_finder_config.io_params.input_params.germline_input,
                                                         config.vj_finder_config.algorithm_params.germline_params);
        auto v_gene_database = db_generator.GenerateVariableDb();
        auto j_gene_database = db_generator.GenerateJoinDb();
        auto v_labeling = cdr_labeler::GermlineDbLabeler(v_gene_database, config.cdrs_params).ComputeLabeling();
        auto j_labeling = cdr_labeler::GermlineDbLabeler(j_gene_database, config.cdrs_params).ComputeLabeling();
        filtered_v_db = v_labeling.CreateFilteredDb();
        filtered_j_db = j_labeling.CreateFilteredDb();
    }
};

std::string GetCDRSequence(const germline_utils::ImmuneGene &immune_gene, annotation_utils::CDRRange cdr_range) {
    auto seqan_cdr = seqan::infixWithLength(immune_gene.seq(), cdr_range.start_pos, cdr_range.length());
    std::stringstream ss;
    ss << seqan_cdr;
    return ss.str();
}

TEST_F(CDRLabelerTest, VGermlineCDRsAreIgBlastConsistent) {
    auto ighv1_18 = filtered_v_db[0];
    INFO("Testing CDR1 and CDR2 for gene " << ighv1_18.name());
    auto v_labeling = cdr_labeler::GermlineDbLabeler(filtered_v_db, config.cdrs_params).ComputeLabeling();
    auto v_cdrs = v_labeling.GetLabelingByGene(ighv1_18);
    std::string computed_cdr1 = GetCDRSequence(ighv1_18, v_cdrs.cdr1);
    std::string computed_cdr2 = GetCDRSequence(ighv1_18, v_cdrs.cdr2);
    std::string igblast_cdr1 = "GGTTACACCTTTACCAGCTATGGT";
    std::string igblast_cdr2 = "ATCAGCGCTTACAATGGTAACACA";
    ASSERT_EQ(igblast_cdr1, computed_cdr1);
    ASSERT_EQ(igblast_cdr2, computed_cdr2);
}

//TEST_F(CDRLabelerTest, JGermlineCDRsAreIgBlastConsistent) {
//    auto ighj = filtered_j_db[0];
//    INFO("Testing CDR3 end for gene " << ighj.name());
//    //auto j_cdrs = j_labeling.GetLabelingByGene(ighj);
//
//}

void CheckFirstAnnotatedClone(const annotation_utils::CDRAnnotatedCloneSet &clone_set) {
    std::string cdr1 = core::seqan_string_to_string(clone_set[0].CDR1());
    ASSERT_EQ(cdr1, std::string("GGATTCACCTTTGATGATTATGGC"));
    std::string cdr2 = core::seqan_string_to_string(clone_set[0].CDR2());
    ASSERT_EQ(cdr2, std::string("ATTAATTGGAATGGTGGTAGCACA"));
    std::string cdr3 = core::seqan_string_to_string(clone_set[0].CDR3());
    ASSERT_EQ(cdr3, std::string("GCGAGAGATCATGATAGTAGTAGCCCGGGGTCCAACTGGTTCGACCCC"));
    std::string aa_seq = core::seqan_string_to_string(clone_set[0].AA());
    ASSERT_EQ(aa_seq, std::string("EVQLVESGGGVVRPGGSLRLSCAASGFTFDDYGMSWVRQAPGKGLEWVSGINWNGGSTGYADSVKGRFTISRDNAKNSLYLQMNSLRAEDTALYHCARDHDSSSPGSNWFDPWGQGTLVTVSS"));
    ASSERT_EQ(clone_set[0].VSHMs().size(), 0);
    ASSERT_EQ(clone_set[0].JSHMs().size(), 0);
}

void CheckSecondAnnotatedClone(const annotation_utils::CDRAnnotatedCloneSet &clone_set) {
    std::string cdr1 = core::seqan_string_to_string(clone_set[1].CDR1());
    ASSERT_EQ(cdr1, std::string("GGATTCACCTTCAGTAGCTATGCT"));
    std::string cdr2 = core::seqan_string_to_string(clone_set[1].CDR2());
    ASSERT_EQ(cdr2, std::string("ATTAATTGGAATGGTGGTAGCACA"));
    std::string cdr3 = core::seqan_string_to_string(clone_set[1].CDR3());
    ASSERT_EQ(cdr3, std::string("GCGAGGAGATTTTCAGAAGGAGCTTTTGATATC"));
    ASSERT_EQ(clone_set[1].VSHMs().size(), 16);
    // ASSERT_EQ(clone_set[1].JSHMs().size(), 2);
    ASSERT_EQ(clone_set[1].JSHMs().size(), 0);
}

TEST_F(CDRLabelerTest, ReadCDRsAreIgBlastConsistent) {
    using namespace cdr_labeler;
    auto v_labeling = GermlineDbLabeler(filtered_v_db, config.cdrs_params).ComputeLabeling();
    INFO("CDR labeling for J gene segments");
    auto j_labeling = GermlineDbLabeler(filtered_j_db, config.cdrs_params).ComputeLabeling();
    INFO("Creation of labeled V and J databases");
    auto labeled_v_db = v_labeling.CreateFilteredDb();
    INFO("Labeled DB of V segments consists of " << labeled_v_db.size() << " records");
    auto labeled_j_db = j_labeling.CreateFilteredDb();
    INFO("Labeled DB of J segments consists of " << labeled_j_db.size() << " records");
    INFO("Alignment against VJ germline segments");
    core::ReadArchive read_archive("test_dataset/cdr_labeler_minitest.fastq");
    vj_finder::VJParallelProcessor processor(read_archive, config.vj_finder_config.algorithm_params,
                                             labeled_v_db, labeled_j_db,
                                             config.run_params.num_threads);
    vj_finder::VJAlignmentInfo alignment_info = processor.Process();
    ReadCDRLabeler read_labeler(config.shm_params, v_labeling, j_labeling);
    auto annotated_clone_set = read_labeler.CreateAnnotatedCloneSet(alignment_info);
    ASSERT_EQ(annotated_clone_set.size(), 2);
    CheckFirstAnnotatedClone(annotated_clone_set);
    CheckSecondAnnotatedClone(annotated_clone_set);
}
