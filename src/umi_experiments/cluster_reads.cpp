#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <segfault_handler.hpp>
#include "utils.hpp"
#include "umi_utils.hpp"
#include "../fast_ig_tools/banded_half_smith_waterman.hpp"

bool readArgs(int argc, char **argv, std::string& reads_file, std::string& output_file) {
    namespace po = boost::program_options;
    po::options_description cmdl_options("Is this needed?");
    cmdl_options.add_options()
            ("help,h", "print help message")
            ("reads,r", po::value<std::string>(&reads_file)->required(), "input file with reads")
            ("output,o", po::value<std::string>(&output_file)->default_value(""), "output file with stats")
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
    Stats(std::vector<size_t> hamming_dist_distribution, std::vector<size_t> sw_dist_distribution)
            : hamming_dist_distribution_(hamming_dist_distribution), sw_dist_distribution_(sw_dist_distribution) {};

    std::string ToString();

    static Stats GetStats(const std::vector<seqan::Dna5String>& input_reads, const std::unordered_map<Umi, std::vector<size_t> >& umi_to_reads);

private:
    std::vector<size_t> hamming_dist_distribution_;
    std::vector<size_t> sw_dist_distribution_;
};

Stats Stats::GetStats(const std::vector<seqan::Dna5String>& input_reads, const std::unordered_map<Umi, std::vector<size_t> >& umi_to_reads) {
    size_t max_read_length = 0;
    for (auto& read : input_reads) {
        max_read_length = std::max(max_read_length, length(read));
    }

    auto get_hamming_dist = [&](const seqan::Dna5String& s1, const seqan::Dna5String& s2) -> size_t {
        return static_cast<size_t>(-half_sw_banded(s1, s2, 0, -1, -1, [](int) -> int { return 0; }, 0));
    };

    std::vector<size_t> hamming_dist_distribution(max_read_length + 1);
    std::vector<size_t> sw_dist_distribution(max_read_length + 1);
    for (auto& entry : umi_to_reads) {
        auto& reads = entry.second;
        for (size_t i = 0; i < reads.size(); i++) {
            for (size_t j = 0; j < i; j++) {
                auto& first = input_reads[reads[i]];
                auto& second = input_reads[reads[j]];
                size_t hamming_dist = get_hamming_dist(first, second);
                size_t sw_dist = get_sw_dist(first, second);
                hamming_dist_distribution[hamming_dist]++;
                sw_dist_distribution[sw_dist]++;
            }
        }
    }
    return Stats(hamming_dist_distribution, sw_dist_distribution);
}

std::string Stats::ToString() {
    size_t max_dist = 0;
    for (size_t i = 0; i < hamming_dist_distribution_.size(); i ++) {
        if (hamming_dist_distribution_[i] != 0 || sw_dist_distribution_[i] != 0) {
            max_dist = i;
        }
    }

    std::stringstream ss;
    char buf[100];
    sprintf(buf, "%10s%20s%20s", "dist", "hamming", "sw");
    ss << buf << std::endl;
    for (size_t i = 0; i <= max_dist; i ++) {
        sprintf(buf, "%10d%20d%20d", static_cast<int>(i), static_cast<int>(hamming_dist_distribution_[i]), static_cast<int>(sw_dist_distribution_[i]));
        ss << buf << std::endl;
    }
    return std::string(buf);
}

int main(int argc, char **argv) {
    segfault_handler sh;
    create_console_logger();
    std::string reads_file;
    std::string output_file;
    if (!readArgs(argc, argv, reads_file, output_file)) {
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
    INFO(umi_to_reads.size() << " unique barcodes found");

    INFO("Calculating stats");
    auto stats = Stats::GetStats(input_reads, umi_to_reads);

    auto string_stats = stats.ToString();
    if (output_file != "") {
        INFO("Saving results");
        std::ofstream output(output_file);
        output << stats.ToString() << std::endl;
        INFO("Stats printed to " << output_file);
    } else {
        std::cout << stats.ToString();
    }
}
