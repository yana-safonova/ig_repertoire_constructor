#include <chrono>
#include <atomic>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/format.hpp>
using bformat = boost::format;

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::CharString;

#include "ig_matcher.hpp"
#include "banded_half_smith_waterman.hpp"
#include "ig_final_alignment.hpp"
#include "fast_ig_tools.hpp"

template<typename T>
size_t find_candidates_num(const T &read,
                           const KmerIndex &kmer2reads,
                           int tau, size_t K) {
    int strategy_int = 1;
    size_t required_read_length = (strategy_int != 0) ? (K * (tau + strategy_int)) : 0;
    if (length(read) < required_read_length) {
        return 0;
    }

    std::vector<int> costs;

    auto hashes = polyhashes(read, K);

    costs.reserve(hashes.size());
    for (size_t hash : hashes) {
        auto it = kmer2reads.find(hash);
        if (it != kmer2reads.cend()) // TODO check it
            costs.push_back(it->second.size());
        else
            costs.push_back(0);
    }

    std::vector<size_t> ind = optimal_coverage(costs, K, tau + 1);

    size_t result = 0;

    for (size_t i : ind) {
        size_t hash = hashes[i];
        auto it = kmer2reads.find(hash);
        if (it != kmer2reads.cend()) {
            result += it->second.size();
        }
    }

    return result;
}

template<typename T>
size_t complexityEstimation(const std::vector<T> &input_reads,
                            const KmerIndex &kmer2reads,
                            int tau,
                            int K) {
    std::atomic<size_t> result;
    result = 0;

    Graph g(input_reads.size());

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 8))
    for (size_t j = 0; j < input_reads.size(); ++j) {
        result += find_candidates_num(input_reads[j], kmer2reads, tau, K);
    }

    return result;
}

bool parse_cmd_line_arguments(int argc, char **argv,
                              std::string &input_file, std::string &output_file,
                              int &tau,
                              int &nthreads) {
    std::string config_file = "";

    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
            ("version,v", "print version string")
            ("help,h", "produce help message")
            ("config,c", po::value<std::string>(&config_file)->default_value(config_file),
             "name of a file of a configuration")
            ("input-file,i", po::value<std::string>(&input_file),
             "name of an input file (FASTA|FASTQ)")
            ("output-file,o", po::value<std::string>(&output_file)->default_value(output_file),
             "file for outputted stats")
            ;

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
            ("tau", po::value<int>(&tau)->default_value(tau),
             "maximum distance value for truncated dist-graph construction")
            ("threads,t", po::value<int>(&nthreads)->default_value(nthreads),
             "the number of parallel threads")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
            ("help-hidden", "show all options, including developers' ones")
            ;

    po::options_description cmdline_options("All command line options");
    cmdline_options.add(generic).add(config).add(hidden);

    po::options_description config_file_options;
    config_file_options.add(config).add(hidden);

    po::options_description visible("Allowed options");
    visible.add(generic).add(config);

    po::positional_options_description p;
    p.add("input-file", -1);
    // p.add("output-file", -1);

    po::variables_map vm;
    store(po::command_line_parser(argc, argv).
            options(cmdline_options).positional(p).run(), vm);
    notify(vm);

    if (config_file != "") {
        std::ifstream ifs(config_file.c_str());
        if (!ifs) {
            cout << "can not open config file: " << config_file << "\n";
            return false;
        } else {
            store(parse_config_file(ifs, config_file_options), vm);
            // reparse cmd line again for update config defaults
            store(po::command_line_parser(argc, argv).
                    options(cmdline_options).positional(p).run(), vm);
            notify(vm);
        }
    }

    if (vm.count("help-hidden")) {
        cout << cmdline_options << std::endl;
        return false;
    }

    if (vm.count("help") || !vm.count("input-file")) { // TODO Process required arguments by the proper way
        cout << visible << "\n";
        return false;
    }

    if (vm.count("version")) {
        cout << "<Some cool name> version 0.1" << vm.count("version") << std::endl;
        return false;
    }

    return true;
}

int main(int argc, char **argv) {
    segfault_handler sh;
    perf_counter pc;
    create_console_logger("");

    INFO("Command line: " << join_cmd_line(argc, argv));

    int tau = 3;
    int nthreads = 4;
    std::string input_file = "cropped.fa";
    std::string output_file = "compl_stats.txt";

    try {
        if (!parse_cmd_line_arguments(argc, argv, input_file, output_file, tau, nthreads)) {
            return 0;
        }
    } catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    INFO("Input reads: " << input_file);

    SeqFileIn seqFileIn_input(input_file.c_str());
    std::vector<CharString> input_ids;
    std::vector<Dna5String> input_reads;

    INFO("Reading input reads starts");
    readRecords(input_ids, input_reads, seqFileIn_input);
    INFO(input_reads.size() << " reads were extracted from " << input_file);

    size_t min_L = 999999999;
    for (const auto read : input_reads) {
        min_L = std::min(min_L, length(read));
    }

    INFO("Minimal length: " << min_L);

    std::ofstream out(output_file);

    for (int K : {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90}) {

        INFO("K-mer index construction. K = " << K);
        auto kmer2reads = kmerIndexConstruction(input_reads, K);

        omp_set_num_threads(nthreads);
        INFO(bformat("Truncated distance graph construction using %d threads starts") % nthreads);
        INFO("Construction of candidates graph");

        size_t complexity = complexityEstimation(input_reads,
                                                 kmer2reads,
                                                 tau, K);

        // Output
        out << K << " " << complexity << " " << static_cast<double>(complexity) / input_reads.size() << std::endl;
    }

    INFO("Stats was written to " << output_file);

    INFO("Running time: " << running_time_format(pc));

    return 0;
}

// vim: ts=4:sw=4
