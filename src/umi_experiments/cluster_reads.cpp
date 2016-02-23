#include <atomic>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <openmp_wrapper.h>
#include <segfault_handler.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include "utils.hpp"
#include "umi_utils.hpp"
#include "../fast_ig_tools/banded_half_smith_waterman.hpp"

bool readArgs(int argc, char **argv, std::string& reads_file, std::string& output_dir, unsigned& thread_count) {
    namespace po = boost::program_options;
    po::options_description cmdl_options("Is this needed?");
    cmdl_options.add_options()
            ("help,h", "print help message")
            ("reads,r", po::value<std::string>(&reads_file)->required(), "input file with reads")
            ("output,o", po::value<std::string>(&output_dir)->default_value(""), "output directory")
            ("threads,t", po::value<unsigned>(&thread_count)->default_value(1), "number of running threads")
            ;
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(cmdl_options).run(), vm);
    if (vm.count("help") || argc == 1) {
        std::cout << cmdl_options << std::endl;
        return false;
    }
    po::notify(vm);
    return true;
}

void group_reads_by_umi(const std::vector<seqan::Dna5String>& umis, std::unordered_map<Umi, std::vector<size_t> >& umi_to_reads) {
    for (size_t i = 0; i < umis.size(); i ++) {
        Umi umi(umis[i]);
        auto& reads = umi_to_reads[umi];
        reads.push_back(i);
    }
}

size_t get_sw_dist(const seqan::Dna5String& first, const seqan::Dna5String& second) {
    vector<size_t> cur(length(first) + 1);
    vector<size_t> prev(length(first) + 1);
    size_t INF = std::numeric_limits<size_t>::max() / 2;
    fill(cur.begin() + 1, cur.end(), INF);
    for (size_t i = 0; i < length(second); i ++) {
        std::swap(cur, prev);
        fill(cur.begin(), cur.end(), INF);
        cur[0] = prev[0] + 1;
        for (size_t j = 1; j <= length(first); j ++) {
            cur[j] = min({
                                 cur[j - 1] + 1,
                                 prev[j] + 1,
                                 prev[j - 1] + (first[j - 1] == second[i] ? 0 : 1)
                         });
        }
    }
    return cur[length(first)];
}

class Stats {
public:
    Stats(std::map<size_t, std::map<size_t, size_t>> hamming_dist_distribution, std::map<size_t, std::map<size_t, size_t>> sw_dist_distribution)
            : hamming_dist_distribution_(hamming_dist_distribution), sw_dist_distribution_(sw_dist_distribution) {};

    std::vector<size_t> GetSizes();

    std::string ToString(size_t umi_size);

    static Stats GetStats(const std::vector<seqan::Dna5String>& input_reads, const std::unordered_map<Umi, std::vector<size_t> >& umi_to_reads, unsigned thread_count);

private:
    std::map<size_t, std::map<size_t, size_t>> hamming_dist_distribution_;
    std::map<size_t, std::map<size_t, size_t>> sw_dist_distribution_;
};

Stats Stats::GetStats(const std::vector<seqan::Dna5String>& input_reads, const std::unordered_map<Umi, std::vector<size_t> >& umi_to_reads, unsigned thread_count) {
    size_t max_read_length = 0;
    for (auto& read : input_reads) {
        max_read_length = std::max(max_read_length, length(read));
    }

    auto get_hamming_dist = [&](const seqan::Dna5String& s1, const seqan::Dna5String& s2) -> size_t {
        return static_cast<size_t>(-half_sw_banded(s1, s2, 0, -1, -1, [](int) -> int { return 0; }, 0));
    };

    std::vector<std::vector<size_t>> read_groups;
    for (auto& entry : umi_to_reads) {
        read_groups.push_back(entry.second);
    }
    std::map<size_t, std::map<size_t, size_t>> hamming_dist_distribution;
    std::map<size_t, std::map<size_t, size_t>> sw_dist_distribution;
    std::atomic<size_t> processed;
    processed = 0;
    size_t next_percent = 1;
    omp_set_num_threads(thread_count);
#pragma omp parallel for schedule(dynamic)
    for (size_t group = 0; group < read_groups.size(); group ++) {
        auto& reads = read_groups[group];
        for (size_t i = 0; i < reads.size(); i++) {
            for (size_t j = 0; j < i; j++) {
                auto& first = input_reads[reads[i]];
                auto& second = input_reads[reads[j]];
                size_t hamming_dist = get_hamming_dist(first, second);
                size_t sw_dist = get_sw_dist(first, second);
#pragma omp critical
                {
                    hamming_dist_distribution[reads.size()][hamming_dist] ++;
                    sw_dist_distribution[reads.size()][sw_dist] ++;
                }
            }
        }
        processed += reads.size();
        while (processed * 100 >= input_reads.size() * next_percent) {
#pragma omp critical
            if (processed * 100 >= input_reads.size() * next_percent) {
                INFO(next_percent << "% of " << input_reads.size() << " reads processed");
                next_percent++;
//                std::cout << Stats(hamming_dist_distribution, sw_dist_distribution).ToString();
            }
        }
    }
    return Stats(hamming_dist_distribution, sw_dist_distribution);
}

std::vector<size_t> Stats::GetSizes() {
    std::vector<size_t> sizes(hamming_dist_distribution_.size());
    size_t current = 0;
    for (auto& entry : hamming_dist_distribution_) {
        sizes[current ++] = entry.first;
    }
    return sizes;
}

std::string Stats::ToString(size_t umi_size) {
    size_t max_dist = 0;
    for (size_t i = 0; i < hamming_dist_distribution_[umi_size].size(); i ++) {
        if (hamming_dist_distribution_[umi_size][i] != 0 || sw_dist_distribution_[umi_size][i] != 0) {
            max_dist = i;
        }
    }

    std::stringstream ss;
    char buf[100];
    sprintf(buf, "%5s%10s%10s", "dist", "hamming", "sw");
    ss << buf << std::endl;
    for (size_t i = 0; i <= max_dist; i ++) {
        sprintf(buf, "%5d%10d%10d", static_cast<int>(i), static_cast<int>(hamming_dist_distribution_[umi_size][i]), static_cast<int>(sw_dist_distribution_[umi_size][i]));
        ss << buf << std::endl;
    }
    return ss.str();
}

void print_umi_reads_distribution_by_size(unordered_map<Umi, vector<size_t>>& umi_to_reads) {
    std::map<size_t, size_t> size_to_count;
    for (auto& entry : umi_to_reads) {
        size_to_count[entry.second.size()] ++;
    }
    for (auto& entry : size_to_count) {
        std::cout << "size " << entry.first << ", count " << entry.second << std::endl;
    }
}

int main(int argc, char **argv) {
    segfault_handler sh;
    create_console_logger();
    std::string reads_file;
    std::string output_dir;
    unsigned thread_count;
    if (!readArgs(argc, argv, reads_file, output_dir, thread_count)) {
        return 0;
    }

    INFO("Reading fastq");
    seqan::SeqFileIn seqFileIn_input(reads_file.c_str());
    std::vector<seqan::CharString> input_ids;
    std::vector<seqan::Dna5String> input_reads;
    readRecords(input_ids, input_reads, seqFileIn_input);
    INFO(input_ids.size() << " records read");

    INFO("Extracting barcodes");
    std::vector<seqan::Dna5String> umis;
    {
        std::vector<seqan::DnaQString> umi_quals;
        extract_barcodes_from_read_ids(input_ids, umis, umi_quals);
    }

    INFO("Grouping reads by barcodes");
    std::unordered_map<Umi, std::vector<size_t> > umi_to_reads;
    group_reads_by_umi(umis, umi_to_reads);
    INFO(umi_to_reads.size() << " unique barcodes found, stats by sizes:");
    print_umi_reads_distribution_by_size(umi_to_reads);

    INFO("Calculating stats");
    auto stats = Stats::GetStats(input_reads, umi_to_reads, thread_count);

    if (output_dir != "") {
        INFO("Saving results");
        if (!boost::filesystem::exists(output_dir)) {
            INFO("Creating directory " << output_dir);
            boost::filesystem::create_directory(output_dir);
        }
        for (auto umi_size : stats.GetSizes()) {
            const auto file_path = boost::filesystem::path(output_dir).append(std::to_string(umi_size));
            std::ofstream output(file_path.string());
            output << stats.ToString(umi_size) << std::endl;
        }
        INFO("Stats printed to " << output_dir);
    } else {
        for (auto umi_size : stats.GetSizes()) {
            std::cout << stats.ToString(umi_size);
        }
    }
}
