#include <vector>
#include <string>
#include <algorithm>
#include <chrono>

#include <build_info.hpp>

#include <iostream>
using std::cout;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "fast_ig_tools.hpp"
#include "ig_trie_compressor.hpp"
#include "utils.hpp"

using fast_ig_tools::Compressor;

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
    std::string idmap_file_name = "";
    bool ignore_tails = true;
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
            ("idmap,m", po::value<std::string>(&idmap_file_name)->default_value(idmap_file_name),
             "map file name; empty (default) for non-producing")
            ;

        // Declare a group of options that will be
        // allowed both on command line and in
        // config file
        po::options_description config("Configuration");
        config.add_options()
            ("ignore-tails,T", po::value<bool>(&ignore_tails)->default_value(ignore_tails),
             "wheather to ignore extra tail of the longest read during read comparison")
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

    INFO("Compression of reads starts");
    auto indices = Compressor::compressed_reads_indices(input_reads,
                                                        ignore_tails ? Compressor::Type::TrieCompressor : Compressor::Type::HashCompressor);
    INFO("Compression of reads finished")

    std::vector<size_t> abundances(indices.size());
    std::vector<size_t> index2newindex(indices.size());

    size_t count = 0;
    for (size_t i = 0; i < indices.size(); ++i) {
        if (indices[i] == i) {
            index2newindex[i] = count;
            ++count;
        }

        abundances[indices[i]] += 1;  // TODO parse input read abundances and add their values here
    }

    for (size_t i = 0; i < indices.size(); ++i) {
        if (indices[i] == i) {
            std::string id = seqan::toCString(input_ids[i]);
            id += "___size___" + std::to_string(abundances[i]);

            seqan::writeRecord(seqFileOut_output, id, input_reads[i]);
        }
    }

    INFO(count << " compressed reads were written to " << output_file);

    if (idmap_file_name != "") {
        std::ofstream idmap_file(idmap_file_name.c_str());

        for (size_t i : indices) {
            idmap_file << index2newindex[i] << "\n";
        }

        INFO("Map from input reads to compressed reads was written to " << idmap_file_name);
    }
    INFO("Running time: " << running_time_format(pc));
    return 0;
}

// vim: ts=4:sw=4
