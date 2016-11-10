#include <utils/io.hpp>
#include <unordered_map>
#include <umi_utils.hpp>
#include <unordered_set>
#include <logger/logger.hpp>
#include "error_analyzer.hpp"

void ErrorAnalyzer::readData(std::string input_file_path) {
    read_seqan_records(input_file_path, ids_, reads_);
}

void ErrorAnalyzer::performAnalysis() {
    std::vector<seqan::Dna5String> barcodes;
    extract_barcodes_from_read_ids(ids_, barcodes);
    std::unordered_map<Umi, std::vector<size_t> > umi_to_reads;
    group_reads_by_umi(barcodes, umi_to_reads);

    for (const auto& entry : umi_to_reads) {
        const auto& barcode = entry.first;
        const auto& read_idx_list = entry.second;
        std::unordered_multiset<std::string> read_to_count;
        for (size_t idx : read_idx_list) {
            read_to_count.insert(seqan_string_to_string(reads_[idx]));
        }
        std::vector<size_t> sizes(read_to_count.size());
        std::transform(read_to_count.begin(), read_to_count.end(), sizes.begin(), [&read_to_count](std::string read) { return read_to_count.count(read); });
        std::sort(sizes.begin(), sizes.end(), std::greater<size_t>());
        VERIFY(!sizes.empty());
        size_t sum = std::accumulate(sizes.begin(), sizes.end(), (size_t) 0);
        VERIFY(read_idx_list.size() == sum);
        if (sizes[0] < 5 || sizes[0] * 10 < sum) continue;
        std::stringstream sstr;
        for (size_t size : sizes) {
            sstr << size << ", ";
        }
        INFO(sstr.str());
    }
}
