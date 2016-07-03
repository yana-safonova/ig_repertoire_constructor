#include <gtest/gtest.h>
#include <logger/log_writers.hpp>
#include <read_archive.hpp>

#include <boost/algorithm/string.hpp>
#include <unordered_set>

#include "../vj_finder/germline_db_generator.hpp"
#include "../vj_finder/vj_alignment_info.hpp"
#include "../vj_finder/vj_query_processing.hpp"

template<typename T>
T string_to_number(std::string str) {
    std::stringstream ss;
    ss << str;
    T res;
    ss >> res;
    return res;
}

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

    AlignmentInfoRecord(std::string line) {
        std::vector<std::string> splits;
        boost::split(splits, line, boost::is_any_of(","));
        read_name = splits[0];
        v_start = string_to_number<size_t>(splits[1]) - 1;
        v_end = string_to_number<size_t>(splits[2]);
        v_score = string_to_number<double>(splits[3]);
        v_gene_name = splits[4];
        j_start = string_to_number<size_t>(splits[5]) - 1;
        j_end = string_to_number<size_t>(splits[6]);
        j_score = string_to_number<double>(splits[7]);
        j_gene_name = splits[8];
    }
};

class AlignmentInfo {
    std::vector<AlignmentInfoRecord> records_;
    std::unordered_map<std::string, size_t> read_names_;
public:
    AlignmentInfo() {}

    void Initialize(std::string csv_filename) {
        path::check_existence(csv_filename);
        std::ifstream fhandler(csv_filename);
        size_t index = 0;
        while(!fhandler.eof()) {
            std::string tmp;
            std::getline(fhandler, tmp);
            if(tmp != "" and index != 0) {
                AlignmentInfoRecord record(tmp);
                records_.push_back(record);
                read_names_[record.read_name] = records_.size() - 1;
            }
            index += 1;
        }
        INFO(records_.size() << " alignment records were extracted from " << csv_filename);
        fhandler.close();
    }

    typedef std::vector<AlignmentInfoRecord>::iterator alignment_info_iterator;

    alignment_info_iterator begin() { return records_.begin(); }

    alignment_info_iterator end() { return records_.end(); }

    size_t size() const { return records_.size(); }

    bool ContainsRead(std::string read_name) const {
        return read_names_.find(read_name) != read_names_.end();
    }

    AlignmentInfoRecord GetRecordByReadName(std::string read_name) const {
        VERIFY_MSG(ContainsRead(read_name), "No alignment info for read " << read_name);
        return records_[read_names_.at(read_name)];
    }
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
double epsilon = 0.01;

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

void CheckMatchingAlignedReadNames(const vj_finder::VJAlignmentInfo &alignment_info) {
    INFO("Check for matching of aligned read names");
    bool read_names_matched = true;
    for(size_t i = 0; i < alignment_info.NumVJHits(); i++)
        if(!old_alignment_info.ContainsRead(alignment_info.GetVJHitsByIndex(i).Read().name)) {
            read_names_matched = false;
            break;
        }
    ASSERT_EQ(read_names_matched, true);
}

bool VGeneNameMatches(AlignmentInfoRecord old_record, vj_finder::VJHits vj_hits) {
    if(old_record.v_gene_name != std::string(seqan::toCString(vj_hits.GetVHitByIndex(0).ImmuneGene().name()))) {
        std::cout << vj_hits.Read().name << std::endl;
        std::cout << old_record.v_gene_name << " - " << vj_hits.GetVHitByIndex(0).ImmuneGene().name() << std::endl;
    }
    return old_record.v_gene_name == std::string(seqan::toCString(vj_hits.GetVHitByIndex(0).ImmuneGene().name()));
}

bool VStartsMatches(AlignmentInfoRecord old_record, vj_finder::VJHits vj_hits) {
    return old_record.v_start == vj_hits.GetVHitByIndex(0).Start();
}

bool VEndsMatches(AlignmentInfoRecord old_record, vj_finder::VJHits vj_hits) {
    return old_record.v_end == vj_hits.GetVHitByIndex(0).End();
}

bool VScoresMatches(AlignmentInfoRecord old_record, vj_finder::VJHits vj_hits) {
    return fabs(old_record.v_score - vj_hits.GetVHitByIndex(0).Score()) <= epsilon;
}

bool JGeneNameMatches(AlignmentInfoRecord old_record, vj_finder::VJHits vj_hits) {
    return old_record.j_gene_name == std::string(seqan::toCString(vj_hits.GetJHitByIndex(0).ImmuneGene().name()));
}

bool JStartsMatches(AlignmentInfoRecord old_record, vj_finder::VJHits vj_hits) {
    return old_record.j_start == vj_hits.GetJHitByIndex(0).Start();
}

bool JEndsMatches(AlignmentInfoRecord old_record, vj_finder::VJHits vj_hits) {
    return old_record.j_end == vj_hits.GetJHitByIndex(0).End();
}

bool JScoresMatches(AlignmentInfoRecord old_record, vj_finder::VJHits vj_hits) {
    return fabs(old_record.j_score - vj_hits.GetJHitByIndex(0).Score()) <= epsilon;
}

bool TwoAlignmentRecordsMatches(AlignmentInfoRecord old_record, vj_finder::VJHits vj_hits) {
    return VGeneNameMatches(old_record, vj_hits) and VStartsMatches(old_record, vj_hits) and
            VEndsMatches(old_record, vj_hits) and JGeneNameMatches(old_record, vj_hits) and
            JStartsMatches(old_record, vj_hits) and JEndsMatches(old_record, vj_hits);
}

void CheckTwoAlignmentRecordsMatches(AlignmentInfoRecord old_record, vj_finder::VJHits vj_hits) {
    ASSERT_EQ(VGeneNameMatches(old_record, vj_hits), true);
    ASSERT_EQ(VStartsMatches(old_record, vj_hits), true);
    ASSERT_EQ(VEndsMatches(old_record, vj_hits), true);
    ASSERT_EQ(VScoresMatches(old_record, vj_hits), true);

    ASSERT_EQ(JGeneNameMatches(old_record, vj_hits), true);
    ASSERT_EQ(JStartsMatches(old_record, vj_hits), true);
    ASSERT_EQ(JEndsMatches(old_record, vj_hits), true);
    ASSERT_EQ(JScoresMatches(old_record, vj_hits), true);
}

void CheckMatchingAlignmentResults(const vj_finder::VJAlignmentInfo &alignment_info) {
    INFO("Check for mathing of alignment results");
    bool alignment_results_mathced = true;
    size_t num_not_matched = 0;
    for(size_t i = 0; i < alignment_info.NumVJHits(); i++) {
        auto vj_hits = alignment_info.GetVJHitsByIndex(i);
        AlignmentInfoRecord old_record = old_alignment_info.GetRecordByReadName(vj_hits.Read().name);
        if(!TwoAlignmentRecordsMatches(old_record, vj_hits))
            num_not_matched += 1;
    }
    std::cout << num_not_matched << std::endl;
}

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
    INFO("Check for number of filtered reads");
    ASSERT_EQ(alignment_info.NumVJHits(), old_alignment_info.size());
    CheckMatchingAlignedReadNames(alignment_info);
    CheckMatchingAlignmentResults(alignment_info);
}