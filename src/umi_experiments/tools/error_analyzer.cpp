#include <utils/io.hpp>
#include <unordered_map>
#include <umi_utils.hpp>
#include <unordered_set>
#include <logger/logger.hpp>
#include "error_analyzer.hpp"

void ErrorAnalyzer::readData(std::string input_file_path) {
    read_seqan_records(input_file_path, ids_, reads_);
}

void ErrorAnalyzer::performAnalysis(std::string& error_stats_dir_path) {
    std::vector<seqan::Dna5String> barcodes;
    extract_barcodes_from_read_ids(ids_, barcodes);
    std::unordered_map<Umi, std::vector<size_t>> umi_to_reads;
    group_reads_by_umi(barcodes, umi_to_reads);

    size_t total_barcodes = 0;
    size_t total_reads = 0;
    size_t total_errors = 0;
    std::vector<double> relative_error_positions;
    std::vector<size_t> error_positions;
    std::vector<size_t> error_positions_reversed;
    std::map<size_t, size_t> errors_in_read_to_cnt;

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
        const auto& super_read_itr = std::find_if(read_to_count.begin(), read_to_count.end(), [&read_to_count, &sizes]( std::string read) { return read_to_count.count( read) == sizes[0]; });
        VERIFY(super_read_itr != read_to_count.end());
        std::string super_read = *super_read_itr;
        for (const auto& read : read_to_count) {
            if (read.length() != super_read.length()) continue;
            total_reads ++;
            size_t errors_in_read = 0;
            for (size_t i = 0; i < read.length(); i ++) {
                if (read[i] != super_read[i]) {
                    errors_in_read ++;
                    relative_error_positions.push_back((double) i / (double) super_read.length());
                    error_positions.push_back(i);
                    error_positions_reversed.push_back(super_read.length() - i - 1);
                }
            }
            total_errors += errors_in_read;
            errors_in_read_to_cnt[errors_in_read] ++;
        }
    }

    if (!error_stats_dir_path.empty()) {
        namespace fs = boost::filesystem;
        const fs::path dir_path = boost::filesystem::path(error_stats_dir_path);
        fs::create_directories(dir_path);

        {
            const std::string relative_errors_path = fs::path(dir_path).append("relative.csv").string();
            INFO("Writing relative error positions to " << relative_errors_path);
            std::ofstream ofs_relative(relative_errors_path);
            for (double p : relative_error_positions) {
                ofs_relative << p << "\n";
            }
        }

        {
            std::sort(error_positions.begin(), error_positions.end());
            const std::string errors_path = fs::path(dir_path).append("error_positions.csv").string();
            INFO("Writing error positions to " << errors_path);
            std::ofstream ofs_absolute(errors_path);
            for (size_t p : error_positions) {
                ofs_absolute << p << "\n";
            }
            std::sort(error_positions_reversed.begin(), error_positions_reversed.end());
            const std::string errors_reversed_path = fs::path(dir_path).append("error_positions_reversed.csv").string();
            INFO("Writing reversed error positions to " << errors_reversed_path);
            std::ofstream ofs_absolute_reversed(errors_reversed_path);
            for (size_t p : error_positions_reversed) {
                ofs_absolute_reversed << p << "\n";
            }
        }

        {
            const std::string distr_path = fs::path(dir_path).append("read_error_distr").string();
            INFO("Writing error count distribution to " << distr_path);
            std::ofstream ofs_distr(distr_path);
            size_t max_error = errors_in_read_to_cnt.rbegin()->first;
            for (size_t e = 0; e <= max_error; e ++) {
                ofs_distr << errors_in_read_to_cnt[e] << "\n";
            }
        }
    }

    INFO("Total reads taken: " << total_reads << " from " << total_barcodes << " barcodes. They have " << total_errors
                               << " errors which is " << ((double) total_errors / (double) total_reads) << " per read in average." );
}
