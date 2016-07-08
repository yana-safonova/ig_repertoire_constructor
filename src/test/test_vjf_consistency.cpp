#include <gtest/gtest.h>
#include <logger/log_writers.hpp>
#include <read_archive.hpp>

#include "../vj_finder/germline_db_generator.hpp"
#include "../vj_finder/vj_alignment_info.hpp"
#include "../vj_finder/vj_query_processing.hpp"

struct AlignmentInfoRecord {
    std::string read_name;
    size_t v_start;
    size_t v_end;
    double v_score;
    std::string v_gene_name;
    size_t j_start;
    size_t j_end;
    double j_score;
    std::string j_gene_name;

    AlignmentInfoRecord(std::string) {
        // split and parse
    }
};

class AlignmentInfo {
    std::vector<AlignmentInfoRecord> records_;
public:
    AlignmentInfo() {}

    void Initialize(std::string csv_filename) {
        path::check_existence(csv_filename);
        std::ifstream fhandler(csv_filename);
        size_t index = 0;
        while(!fhandler.eof()) {
            std::string tmp;
            std::getline(fhandler, tmp);
            if(tmp != "" and index != 0)
                records_.push_back(AlignmentInfoRecord(tmp));
            index += 1;
        }
        INFO(records_.size() << " alignment records were extracted from " << csv_filename);
        fhandler.close();
    }

    typedef std::vector<AlignmentInfoRecord>::iterator alignment_info_iterator;

    alignment_info_iterator begin() { return records_.begin(); }

    alignment_info_iterator end() { return records_.end(); }

    size_t size() const { return records_.size(); }
};

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

void disable_fix_fill_crop(vj_finder::VJFinderConfig& config) {
    config.algorithm_params.fix_crop_fill_params.crop_left = false;
    config.algorithm_params.fix_crop_fill_params.crop_right = false;
    config.algorithm_params.fix_crop_fill_params.fill_left = false;
    config.algorithm_params.fix_crop_fill_params.fill_right = false;
    config.algorithm_params.fix_crop_fill_params.enable_fixing = false;
}

AlignmentInfo old_alignment_info;
vj_finder::VJFinderConfig vj_finder_config;

class VjfConsistencyTest : public ::testing::Test {
public:
    void SetUp() {
        create_console_logger();
        std::string csv_filename = "test_dataset/vj_finder_test/old_vj_finder.csv";
        old_alignment_info.Initialize(csv_filename);
        std::string config_fname = "configs/vj_finder/config.info";
        vj_finder::load(vj_finder_config, config_fname);
        disable_fix_fill_crop(vj_finder_config);
    }
};

TEST_F(VjfConsistencyTest, OldAndNewVjFinderResultsAreConsistent) {
    core::ReadArchive read_archive(vj_finder_config.io_params.input_params.input_reads);
    if(vj_finder_config.io_params.output_params.output_details.fix_spaces)
        read_archive.FixSpacesInHeaders();

    vj_finder::GermlineDbGenerator db_generator(vj_finder_config.io_params.input_params.germline_input,
                                                vj_finder_config.algorithm_params.germline_params);
    INFO("Generation of DB for variable segments...");
    germline_utils::CustomGeneDatabase v_db = db_generator.GenerateVariableDb();
    INFO("Generation of DB for join segments...");
    germline_utils::CustomGeneDatabase j_db = db_generator.GenerateJoinDb();
    vj_finder::VJAlignmentInfo alignment_info;
    for(auto it = read_archive.cbegin(); it != read_archive.cend(); it++) {
        vj_finder::VJQueryProcessor vj_query_processor(vj_finder_config.algorithm_params, v_db, j_db);
        auto vj_hits = vj_query_processor.Process(*it);
        if(!vj_hits) {
            alignment_info.UpdateFilteredReads(*it);
        }
        else {
            alignment_info.UpdateHits(*vj_hits);
        }
    }
    ASSERT_EQ(alignment_info.NumVJHits(), old_alignment_info.size());
}