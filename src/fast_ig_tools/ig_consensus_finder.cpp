#include <fstream>
#include <vector>
#include <cassert>
#include <algorithm>

#include <unordered_map>

#include <iostream>
using std::cout;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/format.hpp>
using bformat = boost::format;

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::SeqFileOut;
using seqan::CharString;

#include "ig_final_alignment.hpp"
#include "fast_ig_tools.hpp"


std::pair<std::unordered_map<std::string, size_t>, std::vector<std::string>> read_rcm_file_string(const std::string &file_name) {
    std::ifstream rcm(file_name.c_str());

    std::unordered_map<std::string, size_t> result;
    std::unordered_map<std::string, size_t> barcode2num;

    std::string id, target;
    std::vector<std::string> targets;
    size_t num = 0;
    while (rcm >> id >> target) {
        if (barcode2num.count(target) == 0) { // TODO Use iterator here??
            barcode2num[target] = num++;
            targets.push_back(target);
            assert(targets.size() == num);
        }

        result[id] = barcode2num[target];
    }

    return { result, targets };
}


int main(int argc, char **argv) {
    segfault_handler sh;
    perf_counter pc;
    create_console_logger("");

    INFO("Command line: " << join_cmd_line(argc, argv));

    int nthreads = 4;
    std::string reads_file = "input.fa";
    std::string output_file = "repertoire.fa";
    std::string rcm_file = "cropped.rcm";
    bool use_hamming_alignment = false;
    std::string config_file = "";

    // Parse cmd-line arguments
    try {
        // Declare a group of options that will be
        // allowed only on command line
        po::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version string")
            ("help,h", "produce help message")
            ("config,c", po::value<std::string>(&config_file)->default_value(config_file),
             "name of a file of a configuration")
            ("input-file,i", po::value<std::string>(&reads_file),
             "name of the input file (FASTA|FASTQ)")
            ("output-file,o", po::value<std::string>(&output_file)->default_value(output_file),
             "output file for final repertoire")
            ("rcm-file,R", po::value<std::string>(&rcm_file)->default_value(rcm_file),
             "input RCM-file")
            ("hamming,H",
             "use Hamming-based (position-wise) multiple alignment instead of seqan's one")
            ;

        // Declare a group of options that will be
        // allowed both on command line and in
        // config file
        po::options_description config("Configuration");
        config.add_options()
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


        po::variables_map vm;
        store(po::command_line_parser(argc, argv).
              options(cmdline_options).run(), vm);
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
                      options(cmdline_options).run(), vm);
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

        if (vm.count("hamming")) {
            use_hamming_alignment = true;
        }
    } catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    INFO("Input files: " << reads_file << ", " << rcm_file);

    std::vector<Dna5String> input_reads;
    std::vector<CharString> input_ids;
    SeqFileIn seqFileIn_input(reads_file.c_str());

    INFO("Reading input reads starts");
    readRecords(input_ids, input_reads, seqFileIn_input);
    INFO(input_reads.size() << " reads were extracted from " << reads_file);

    std::vector<size_t> component_indices;
    component_indices.resize(input_reads.size());

    INFO("Reading read-cluster map starts");
    auto rcm = read_rcm_file_string(rcm_file);
    for (size_t i = 0; i < input_reads.size(); ++i) {
        const char *id = toCString(input_ids[i]);
        assert(rcm.first.count(id));
        component_indices[i] = rcm.first[id];
    }

    INFO(rcm.second.size() << " clusters were extracted from " << rcm_file);

    std::vector<std::vector<size_t>> component2id(rcm.second.size());
    for (size_t i = 0; i < component_indices.size(); ++i) {
        component2id[component_indices[i]].push_back(i);
    }

    size_t max_component_size = 0;
    for (const auto &_ : component2id) {
        max_component_size = std::max(max_component_size, _.size());
    }
    INFO(bformat("Size of maximal cluster: %d") % max_component_size);

    auto abundances = find_abundances(input_ids);

    std::vector<Dna5String> consensuses(component2id.size());
    std::vector<size_t> comp_abundances(component2id.size());

    omp_set_num_threads(nthreads);
    INFO(bformat("Computation of consensus using %d threads starts") % nthreads);

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 8))
    for (size_t comp = 0; comp < component2id.size(); ++comp) {
        if (component2id[comp].empty()) {
            continue;
        }

        if (use_hamming_alignment) {
            consensuses[comp] = consensus_hamming_limited_coverage(input_reads, component2id[comp], 1005000);
        } else {
            consensuses[comp] = consensus(input_reads, component2id[comp]);
        }

        for (size_t i : component2id[comp]) {
            comp_abundances[comp] += abundances[i];
        }
    }

    INFO("Saving results");

    SeqFileOut seqFileOut_output(output_file.c_str());

    for (size_t comp = 0; comp < component2id.size(); ++comp) {
        if (component2id[comp].empty()) {
            continue;
        }

        size_t abundance = comp_abundances[comp];

        bformat fmt("cluster___%s___size___%d");
        fmt % rcm.second[comp] % abundance;
        std::string id = fmt.str();

        seqan::writeRecord(seqFileOut_output, id, consensuses[comp]);
    }

    INFO("Final repertoire was written to " << output_file);

    INFO("Running time: " << running_time_format(pc));

    return 0;
}

// vim: ts=4:sw=4
