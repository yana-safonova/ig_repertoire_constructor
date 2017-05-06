#include "fast_ig_tools.hpp"
#include "ig_final_alignment.hpp"
#include "utils.hpp"
#include <boost/program_options.hpp>
#include <build_info.hpp>
#include <iostream>

#include "hogc.hpp"
using hogc::hammingGraph;
using hogc::discardedReads;

#include <seqan/seq_io.h>
#undef NDEBUG

using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::CharString;

struct HGCParameters {
    size_t k = 10;
    size_t tau = 4;
    size_t threads = 4;
    std::string input_file = "";
    std::string output_file = "";
    size_t strategy = 1;
    bool export_abundances = false;
    bool ignore_tails = true;
};

bool parse_cmd_line_arguments(int argc, char **argv, HGCParameters &args) {
    using std::cout;
    using std::cerr;
    using std::endl;
    namespace po = boost::program_options;

    std::string config_file;

    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()("version,v", "print version string")
            ("help,h",
            "produce help message")
            ("config,c", po::value<std::string>(&config_file),
            "name of a file of a configuration")
            ("input-file,i", po::value<std::string>(&args.input_file)->required(),
            "name of an input file (FASTA|FASTQ)")
            ("output-file,o", po::value<std::string>(&args.output_file),
            "file for outputted truncated dist-graph in METIS format")
            ;

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
            ("word-size,k", po::value<size_t>(&args.k)->default_value(args.k),
            "word size for k-mer index construction")
            ("strategy,S", po::value<size_t>(&args.strategy)->default_value(args.strategy),
            "strategy type (0 --- naive, 1 --- single, 2 --- pair, 3 --- triple, etc)")
            ("tau", po::value<size_t>(&args.tau)->default_value(args.tau),
            "maximum distance value for truncated dist-graph construction")
            ("ignore-tails,T", po::value<bool>(&args.ignore_tails)->default_value(args.ignore_tails)->implicit_value(args.ignore_tails),
            "wheather to ignore extra tail of the longest read during read comparison")
            ("threads,t", po::value<size_t>(&args.threads)->default_value(args.threads),
            "the number of parallel threads")
            ("export-abundances,A", "export read abundances to output graph file")
            ("no-export-abundances", "don't export read abundances to output graph file (default)")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()("help-hidden", "show all options, including developers' ones");

    po::options_description cmdline_options("All command line options");
    cmdline_options.add(generic).add(config).add(hidden);

    po::options_description config_file_options;
    config_file_options.add(config).add(hidden);

    po::options_description visible("Allowed options");
    visible.add(generic).add(config);

    po::positional_options_description p;
    p.add("input-file", 1);
    p.add("output-file", 1);

    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

    if (vm.count("help-hidden")) {
        cout << cmdline_options << endl;
        return false;
    }

    if (vm.count("help")) {
        cout << visible << endl;
        return false;
    }

    if (vm.count("version")) {
        cout << bformat("Hamming Graph Constructor, part of IgReC version %s; git version: %s") % build_info::version % build_info::git_hash7 << endl;
        return false;
    }

    if (config_file != "") {
        std::ifstream ifs(config_file.c_str());
        if (!ifs) {
            cout << "can not open config file: " << config_file << "\n";
            return 0;
        } else {
            store(parse_config_file(ifs, config_file_options), vm);
            // reparse cmd line again for update config defaults
            store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);
        }
    }

    try {
        notify(vm);
    } catch (po::error &e) {
        cout << "Parser error: " << e.what() << endl;
    }

    if (vm.count("export-abundances")) {
        args.export_abundances = true;
    }

    if (vm.count("no-export-abundances")) {
        args.export_abundances = false;
    }

    return true;
}

int main(int argc, char **argv) {
    segfault_handler sh;
    perf_counter pc;
    create_console_logger("");

    HGCParameters args;
    try {
        if (!parse_cmd_line_arguments(argc, argv, args)) {
            return 0;
        }
    } catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    INFO("Command line: " << join_cmd_line(argc, argv));
    INFO("Input reads: " << args.input_file);
    INFO("k = " << args.k << ", tau = " << args.tau);

    SeqFileIn seqFileIn_input(args.input_file.c_str());
    std::vector<CharString> input_ids;
    std::vector<Dna5String> input_reads;

    INFO("Reading input reads starts");
    readRecords(input_ids, input_reads, seqFileIn_input);
    INFO(input_reads.size() << " reads were extracted from " << args.input_file);

    std::vector<size_t> abundances;
    if (args.export_abundances) {
        abundances = find_abundances(input_ids);
    }
    // Free memory
    input_ids.clear();

    omp_set_num_threads(args.threads);
    INFO(bformat("Truncated distance graph construction using %d threads starts") % args.threads);
    INFO("Strategy " << args.strategy << " was chosen");

    size_t discarded_reads = discardedReads(input_reads, args.tau, args.k, args.strategy);
    if (discarded_reads) {
        WARN("Reads discarded: "<< discarded_reads);
    }

    size_t num_of_dist_computations;
    auto graph = hammingGraph(input_reads,
                              args.tau,
                              args.ignore_tails,
                              args.k,
                              args.strategy,
                              &num_of_dist_computations);
    INFO("Check graph");
    assert(graph.verify());

    INFO("Simularity computations: " << num_of_dist_computations << ", average " << static_cast<double>(num_of_dist_computations) / input_reads.size() << " per read");

    size_t num_of_edges = graph.nEdges();
    INFO("Edges found: " << num_of_edges);
    INFO("Strategy efficiency: " << static_cast<double>(num_of_edges) / num_of_dist_computations);

    if (args.export_abundances) {
        INFO("Saving graph (with abundances)");
        graph.setVerticesWeights(std::move(abundances));
        graph.writeMetis(args.output_file);
    } else {
        INFO("Saving graph (without abundances)");
        graph.writeMetis(args.output_file);
    }

    INFO("Graph was written to " << args.output_file);

    INFO("Running time: " << running_time_format(pc));

    return 0;
}
