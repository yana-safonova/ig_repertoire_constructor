#include <chrono>
#include <atomic>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#include "fast_ig_tools.hpp"
#include "ig_final_alignment.hpp"
#include "ig_matcher.hpp"
#include "banded_half_smith_waterman.hpp"
#include "utils.hpp"

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::CharString;

template<typename T>
std::pair<size_t, std::vector<size_t>> find_candidates_num(const T &read,
                                                           const KmerIndex &kmer2reads,
                                                           int tau, size_t K) {
    unsigned strategy = 1;
    size_t required_read_length = (strategy != 0) ? (K * (tau + strategy)) : 0;
    if (length(read) < required_read_length) {
        return { 0, {} };
    }

    std::vector<size_t> costs;

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

    return { result, ind };
}

template<typename T>
size_t complexityEstimation(const std::vector<T> &input_reads,
                            const KmerIndex &kmer2reads,
                            int tau,
                            int K,
                            std::vector<std::vector<size_t>> &opt_kmers) {

    Graph g(input_reads.size());

    size_t result = 0;
    SEQAN_OMP_PRAGMA(parallel for reduction(+:result) schedule(dynamic, 8))
    for (size_t j = 0; j < input_reads.size(); ++j) {
        auto _ = find_candidates_num(input_reads[j], kmer2reads, tau, K);
        opt_kmers[j] = _.second;
        result += _.first;
    }

    return result;
}

bool parse_cmd_line_arguments(int argc, char **argv,
                              std::string &input_file, std::string &output_file,
                              int &tau,
                              int &nthreads,
                              bool &save_opt_kmers,
                              int &k_step) {
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
            ("optimal-kmers,O", "save optimal k-mer positions")
            ;

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
            ("tau", po::value<int>(&tau)->default_value(tau),
             "maximum distance value for truncated dist-graph construction")
            ("k-step", po::value<int>(&k_step)->default_value(k_step),
             "step for optimal k-mer size search")
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

    if (vm.count("optimal-kmers")) {
        save_opt_kmers = true;
    } else {
        save_opt_kmers = false;
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
    int k_step = 5;
    int nthreads = 4;
    std::string input_file = "cropped.fa";
    std::string output_file = "compl_stats.txt";
    bool save_opt_kmers = false;

    try {
        if (!parse_cmd_line_arguments(argc, argv, input_file, output_file, tau, nthreads, save_opt_kmers, k_step)) {
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
    for (const auto &read : input_reads) {
        min_L = std::min(min_L, length(read));
    }

    INFO("Minimal length: " << min_L);

    std::ofstream out(output_file);
    out << "# Input file: " << input_file << std::endl;
    out << "# Reads: " << input_reads.size() << std::endl;
    auto bfc = (long long)(input_reads.size()) * (input_reads.size() - 1) / 2;
    out << "# Brute-force d_count: " << bfc << std::endl;
    out << "# Brute-force av_d_count: " << double(input_reads.size() - 1) / 2  << std::endl;
    out << "# Minimal lenght: " << min_L << std::endl;
    out << "# tau: " << tau << std::endl;
    out << "k\td_count\tav_d_count" << std::endl;

    for (int K = 5; K <= std::max(static_cast<int>(min_L) / (tau + 1), 100); K += k_step) {
        INFO("K-mer index construction. K = " << K);
        auto kmer2reads = kmerIndexConstruction(input_reads, K);

        omp_set_num_threads(nthreads);

        std::vector<std::vector<size_t>> opt_kmers(input_reads.size());

        INFO(bformat("Complexity estimation using %d threads starts") % nthreads);
        size_t complexity = complexityEstimation(input_reads,
                                                 kmer2reads,
                                                 tau, K,
                                                 opt_kmers);

        // Optimal k-mers output
        if (save_opt_kmers) {
            bformat fmt("opt_kmers_tau_%d_k_%d.txt");
            fmt % tau % K;
            std::ofstream opt_kmers_out(fmt.str());
            for (size_t j = 0; j < input_reads.size(); ++j) {
                opt_kmers_out << j + 1 << " " << length(input_reads[j]);
                for (const auto &k : opt_kmers[j]) {
                    opt_kmers_out << " " << k + 1;
                }
                opt_kmers_out << std::endl;
            }
        }
        out << K << "\t" << complexity << "\t" << static_cast<double>(complexity) / static_cast<double>(input_reads.size()) << std::endl;
    }

    INFO("Stats was written to " << output_file);

    INFO("Running time: " << running_time_format(pc));

    return 0;
}

// vim: ts=4:sw=4
