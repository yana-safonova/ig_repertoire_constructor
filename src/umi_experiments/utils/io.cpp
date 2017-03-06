#include <unordered_set>
#include <string>

#include <seqan/seq_io.h>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include <logger/logger.hpp>
#include <verify.hpp>
#include "../../ig_tools/utils/string_tools.hpp"
#include <umi_utils.hpp>

bool create_new_directory(boost::filesystem::path path) {
    namespace fs = boost::filesystem;
    bool exists = fs::exists(path);
    if (exists) {
        fs::remove_all(path);
    }
    fs::create_directory(path);
    return !exists;
}

void write_seqan_records(boost::filesystem::path path, std::vector<seqan::CharString> ids, std::vector<seqan::Dna5String> reads) {
    namespace fs = boost::filesystem;
    seqan::SeqFileOut clusters_file(path.string().c_str());
    writeRecords(clusters_file, ids, reads);
}

void read_seqan_records(const std::string& input_file_name, std::vector<seqan::CharString>& ids, std::vector<seqan::Dna5String>& reads) {
    INFO("Reading sequences from " << input_file_name);
    seqan::SeqFileIn reads_file(input_file_name.c_str());
    readRecords(ids, reads, reads_file);
    INFO(reads.size() << " sequences read");
}

std::unordered_map<seqan::CharString, size_t> read_rcm_file(std::string file_path) {
    INFO("Reading rcm from " << file_path);
    std::unordered_map<seqan::CharString, size_t> rcm;
    std::unordered_set<size_t> clusters;
    std::ifstream file(file_path);
    std::string s;
    while (std::getline(file, s)) {
        const auto& mapping = split(s, '\t');
        VERIFY(mapping.size() >= 1 && mapping.size() <= 2);
        if (mapping.size() == 2) {
            size_t cluster = std::stoull(mapping[1]);
            rcm[seqan::CharString(mapping[0])] = cluster;
            clusters.insert(cluster);
        }
    }
    INFO(rcm.size() << " nonempty mappings read. " << clusters.size() << " clusters in total.");
    return rcm;
}
