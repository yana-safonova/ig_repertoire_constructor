#include <utils/io.hpp>
#include <unordered_map>
#include <umi_utils.hpp>
#include <unordered_set>
#include <logger/logger.hpp>
#include "error_analyzer.hpp"

void ErrorAnalyzer::readData(std::string input_file_path) {
    read_seqan_records(input_file_path, ids_, reads_);
}

void ErrorAnalyzer::performAnalysis(std::string& error_pos_file_path) {
    std::vector<seqan::Dna5String> barcodes;
    extract_barcodes_from_read_ids(ids_, barcodes);
    std::unordered_map<Umi, std::vector<size_t> > umi_to_reads;
    group_reads_by_umi(barcodes, umi_to_reads);

    size_t total_barcodes = 0;
    size_t total_reads = 0;
    size_t total_errors = 0;
    std::vector<double> error_positions;

    for (const auto& entry : umi_to_reads) {
        const auto& read_idx_list = entry.second;
        std::unordered_multiset<std::string> read_to_count;
        for (size_t idx : read_idx_list) {
            read_to_count.insert(seqan_string_to_string(reads_[idx]));
        }
        std::unordered_set<std::string> read_set(read_to_count.begin(), read_to_count.end());
        std::vector<size_t> sizes(read_set.size());
        std::transform(read_set.begin(), read_set.end(), sizes.begin(), [&read_to_count](std::string read) { return read_to_count.count(read); });
        std::sort(sizes.begin(), sizes.end(), std::greater<size_t>());
        VERIFY_MSG(!sizes.empty(), "Empty group unexpected.");
        size_t sum = std::accumulate(sizes.begin(), sizes.end(), (size_t) 0);
        VERIFY_MSG(read_idx_list.size() == sum, boost::format("Read group sizes sum to %1% while there were %2% reads sharing barcode.") % sum % read_idx_list.size());
        if (sizes[0] < 5 || sizes[0] * 10 < sum || (sizes.size() >= 2 && sizes[1] * 3 > sizes[0] * 2)) continue;
//        std::stringstream sstr;
//        for (size_t size : sizes) {
//            sstr << size << ", ";
//        }
//        INFO(sstr.str());

        total_barcodes ++;
        std::string super_read = *std::find_if(read_to_count.begin(), read_to_count.end(), [&read_to_count, &sizes](std::string read) { return read_to_count.count(read) == sizes[0]; });
        for (const auto& read : read_to_count) {
            if (read.length() != super_read.length()) continue;
            total_reads ++;
            for (size_t i = 0; i < read.length(); i ++) {
                if (read[i] != super_read[i]) {
                    total_errors ++;
                    error_positions.push_back((double) i / (double) super_read.length());
                }
            }
        }
    }

    if (!error_pos_file_path.empty()) {
        INFO("Writing error positions to " << error_pos_file_path);
        std::ofstream ofs(error_pos_file_path);
        for (double p : error_positions) {
            ofs << p << "\n";
        }
    }

    INFO("Total reads taken: " << total_reads << " from " << total_barcodes << " barcodes. They have " << total_errors
                               << " errors which is " << ((double) total_errors / (double) total_reads) << " per read in average." );
}
