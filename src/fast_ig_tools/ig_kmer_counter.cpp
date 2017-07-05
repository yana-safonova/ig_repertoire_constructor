#include <chrono>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <sstream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#include "fast_ig_tools.hpp"
#include "ig_matcher.hpp"
#include "utils.hpp"

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::CharString;

bool parse_cmd_line_arguments(int argc, char **argv,
                              std::string &input_file,
                              std::string &output_file,
                              int &K) {
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
             "file for outputted k-mer statistics")
            ;

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
            ("word-size,k", po::value<int>(&K)->default_value(K),
             "word size for k-mer index construction")
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

    int K = 36; // anchor length
    std::string input_file = "cropped.fa";
    std::string output_file = "k_mer_stats.txt";

    try {
        if (!parse_cmd_line_arguments(argc, argv, input_file, output_file, K)) {
            return 0;
        }
    } catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    INFO("Input reads: " << input_file);

    INFO("K = " << K);

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

    INFO(input_reads.size() << " reads were extracted from " << input_file);
    INFO("Minimal length: " << min_L);

    INFO("K-mer index construction");
    // auto kmer2reads = kmerIndexConstruction(input_reads, K);


    std::unordered_map<std::string, size_t> prefix_count;
    for (const auto read : input_reads) {
        if (length(read) >= static_cast<size_t>(K)) {
            std::stringstream ss;
            ss << seqan::prefix(read, K);
            std::string s;
            ss >> s;
            prefix_count[s] += 1;
        }
    }
    const auto &kmer2reads = prefix_count;

    std::vector<size_t> kmer_abundances;
    kmer_abundances.reserve(kmer2reads.size());
    for (const auto _ : kmer2reads) {
        // kmer_abundances.push_back(_.second.size());
        kmer_abundances.push_back(_.second);
    }
    std::sort(kmer_abundances.rbegin(), kmer_abundances.rend());

    uint64_t complexity = 0;
    for (size_t abundancy : kmer_abundances) {
        complexity += uint64_t(abundancy) * (abundancy - 1) / 2;
    }

    std::ofstream out(output_file);

    for (size_t abundancy : kmer_abundances) {
        out << abundancy << std::endl;
    }

    INFO("Complexity of naive k-mer-index " << complexity << " dist. computations");
    INFO("Stats was written to " << output_file);
    INFO("Running time: " << running_time_format(pc));

    return 0;
}

// vim: ts=4:sw=4
