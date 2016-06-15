#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
#include <chrono>

#include <build_info.hpp>

#include <iostream>
using std::cout;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "fast_ig_tools.hpp"
#include "ig_trie_compressor.hpp"

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::SeqFileOut;
using seqan::CharString;
using seqan::length;

int main(int argc, char **argv) {
    segfault_handler sh;
    perf_counter pc;
    create_console_logger("");

    std::string input_file = "input.fa";
    std::string output_file = "output.fa";
    std::string rcm_file_name = "";
    try {
        // Declare a group of options that will be
        // allowed only on command line
        po::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version string")
            ("help,h", "produce help message")
            ("config-file,c", "name of a file of a configuration")
            ("input-file,i", po::value<std::string>(&input_file)->required(),
             "name of the input file (FASTA|FASTQ)")
            ("output-file,o", po::value<std::string>(&output_file)->default_value(output_file),
             "name of the output file (FASTA|FASTQ)")
            ("rcm,m", po::value<std::string>(&rcm_file_name)->default_value(rcm_file_name),
             "file name for read-cluster map (rcm); leave empty to skip rcm file generation")
            ;

        // Declare a group of options that will be
        // allowed both on command line and in
        // config file
        po::options_description config("Configuration");
        config.add_options()
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
        p.add("input-file", 1);

        po::variables_map vm;
        store(po::command_line_parser(argc, argv).
              options(cmdline_options).positional(p).run(), vm);


        if (vm.count("help-hidden")) {
            cout << cmdline_options << std::endl;
            return 0;
        }

        if (vm.count("help")) {
            cout << visible << std::endl;
            return 0;
        }

        if (vm.count("version")) {
            cout << bformat("IG Trie compressor, part of IgReC version %s; git version: %s") % build_info::version % build_info::git_hash7 << std::endl;
            return 0;
        }

        if (vm.count("config-file")) {
            std::string config_file = vm["config-file"].as<std::string>();

            std::ifstream ifs(config_file.c_str());
            if (!ifs) {
                cout << "can not open config file: " << config_file << "\n";
                return 0;
            } else {
                store(parse_config_file(ifs, config_file_options), vm);
                // reparse cmd line again for update config defaults
                store(po::command_line_parser(argc, argv).
                      options(cmdline_options).positional(p).run(), vm);
            }
        }

        try {
            notify(vm);
        } catch (po::error &e) {
            cout << "Parser error: " << e.what() << std::endl;
        }
    } catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    INFO("Command line: " << join_cmd_line(argc, argv));
    INFO("Input reads: " << input_file);
    INFO("Output filename: " << output_file);

    SeqFileIn seqFileIn_input(input_file.c_str());
    SeqFileOut seqFileOut_output(output_file.c_str());
    std::vector<CharString> input_ids;
    std::vector<Dna5String> input_reads;

    INFO("Reading input reads starts");
    readRecords(input_ids, input_reads, seqFileIn_input);
    INFO(length(input_reads) << " reads were extracted from " << input_file);

    INFO("Construction of trie starts");
    Trie<seqan::Dna5> trie(input_reads);

    auto representative2abundance = trie.checkout(length(input_reads));
    std::vector<std::pair<size_t, size_t>> representatives(representative2abundance.cbegin(), representative2abundance.cend());
    std::sort(representatives.begin(), representatives.end());
    INFO("Unique prefixes were collected");

    size_t count = 0;
    for (const auto &it : representatives) {
        size_t index = it.first;
        size_t abundance = it.second;

        std::string id = seqan::toCString(input_ids[index]);
        id += "___size___" + std::to_string(abundance);

        seqan::writeRecord(seqFileOut_output, id, input_reads[index]);

        count += abundance;
    }

    assert(count == length(input_reads));

    INFO(representatives.size() << " compressed reads were written to " << output_file);

    if (rcm_file_name != "") {
        std::ofstream rcm_file(rcm_file_name.c_str());
        std::unordered_map<size_t, size_t> representative2compressed_position;
        for (size_t group = 0; group < representatives.size(); group ++) {
            representative2compressed_position[representatives[group].first] = group;
        }
        auto representative2members = trie.checkout_ids(length(input_reads));
        std::vector<size_t> read2compressed_position(length(input_reads));
        for (const auto &entry : representative2members) {
            const size_t representative = entry.first;
            const auto &members = entry.second;
            for (size_t member : members) {
                read2compressed_position[member] = representative2compressed_position[representative];
            }
        }

        for (size_t read = 0; read < input_reads.size(); read ++) {
            rcm_file << input_ids[read] << "\t" << read2compressed_position[read] << "\n";
        }

        INFO("Map from input reads to compressed reads was written to " << rcm_file_name);
    }
    INFO("Construction of trie finished")
    INFO("Running time: " << running_time_format(pc));
    return 0;
}

// vim: ts=4:sw=4
