#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
#include <chrono>

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::SeqFileOut;
using seqan::CharString;
using seqan::length;


#include <boost/format.hpp>
using bformat = boost::format;

#include <iostream>
using std::cout;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "ig_trie_compressor.hpp"
#include "fast_ig_tools.hpp"


int main(int argc, char **argv) {
    segfault_handler sh;
    perf_counter pc;
    create_console_logger("");

    INFO("Command line: " << join_cmd_line(argc, argv));

    std::string input_file = "input.fa";
    std::string output_file = "output.fa";
    std::string idmap_file_name = "";
    try {
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

        po::variables_map vm;
        store(po::command_line_parser(argc, argv).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);

        if (config_file != "") {
            std::ifstream ifs(config_file.c_str());
            if (!ifs) {
                cout << "can not open config file: " << config_file << "\n";
                return 0;
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
            return 0;
        }

        if (vm.count("help") || !vm.count("input-file")) { // TODO Process required arguments by the proper way
            cout << visible << "\n";
            return 0;
        }

        if (vm.count("version")) {
            cout << "<Some cool name> version 0.1" << vm.count("version") << std::endl;
            return 0;
        }
    } catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

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

    auto result__ = trie.checkout(length(input_reads));
    std::vector<std::pair<size_t, size_t>> result(result__.cbegin(), result__.cend());
    std::sort(result.begin(), result.end());
    INFO("Unique prefixes were collected");

    size_t count = 0;
    for (const auto &_ : result) {
        size_t index = _.first;
        size_t abundance = _.second;

        count += abundance;

        std::string id = seqan::toCString(input_ids[index]);
        id += "_abundance:" + std::to_string(abundance);

        seqan::writeRecord(seqFileOut_output, id, input_reads[index]);
    }

    assert(count == length(input_reads));

    INFO(result.size() << " compressed reads were written to " << output_file);

    if (idmap_file_name != "") {
        std::ofstream idmap_file(idmap_file_name.c_str());
        std::vector<size_t> idmap(length(input_reads));
        auto result = trie.checkout_ids(length(input_reads));
        for (const auto &_ : result) {
            size_t index = _.first;
            const auto &ids = _.second;

            for (size_t id : ids) {
                idmap[id] = index;
            }
        }

        // TODO Refactor it ASAP
        std::vector<size_t> ord(length(input_reads), 0);
        for (size_t id : idmap) {
            ord[id] = 1;
        }

        size_t ii = 0;
        for (size_t &o : ord) {
            if (o) {
                o = ii;
                ++ii;
            }
        }


        for (size_t id : idmap) {
            idmap_file << ord[id] << "\n";
        }

        INFO("Map from input reads to compressed reads was written to " << idmap_file_name);
    }
    INFO("Construction of trie finished")
    INFO("Running time: " << running_time_format(pc));
    return 0;
}

// vim: ts=4:sw=4
