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

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::CharString;

#include "ig_matcher.hpp"
#include "banded_half_smith_waterman.hpp"
#include "ig_final_alignment.hpp"
#include "utils.hpp"
#include <build_info.hpp>

struct SWGCParam {
    unsigned k = 10;
    unsigned tau = 4;
    unsigned nthreads = 4;
    std::string input_file = "";
    std::string output_file = "output.graph";
    std::string reference_file = "";
    unsigned strategy = 3;
    unsigned max_indels = 0;
    bool export_abundances = false;
    bool ignore_tails = true;
};


bool parse_cmd_line_arguments(int argc, char **argv, SWGCParam &args) {
    std::string config_file = "";

    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
            ("version,v", "print version string")
            ("help,h", "produce help message")
            ("config,c", po::value<std::string>(&config_file),
             "name of a file of a configuration")
            ("input-file,i", po::value<std::string>(&args.input_file),
             "name of an input file (FASTA|FASTQ)")
            ("reference-file,r", po::value<std::string>(&args.reference_file)->default_value(args.reference_file),
             "name of an input file (FASTA|FASTQ)")
            ("output-file,o", po::value<std::string>(&args.output_file),
             "file for outputted truncated dist-graph in METIS format")
            ("export-abundances,A", "export read abundances to output graph file")
            ("no-export-abundances", "don't export read abundances to output graph file (default)")
            ;

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
            ("word-size,k", po::value<unsigned>(&args.k)->default_value(args.k),
             "word size for k-mer index construction")
            ("ignore-tails,T", po::value<bool>(&args.ignore_tails)->default_value(args.ignore_tails),
             "wheather to ignore extra tail of the longest read during read comparison")
            ("strategy,S", po::value<unsigned>(&args.strategy)->default_value(args.strategy),
             "strategy type (0 --- naive, 1 --- single, 2 --- pair, 3 --- triple, etc)")
            ("tau", po::value<unsigned>(&args.tau)->default_value(args.tau),
             "maximum distance value for truncated dist-graph construction")
            ("max-indels", po::value<unsigned>(&args.max_indels)->default_value(args.max_indels),
             "maximum number of indels in Levenshtein distance")
            ("threads,t", po::value<unsigned>(&args.nthreads)->default_value(args.nthreads),
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
    p.add("input-file", 1);
    p.add("output-file", 1);

    po::variables_map vm;
    store(po::command_line_parser(argc, argv).
          options(cmdline_options).positional(p).run(), vm);


    if (vm.count("help-hidden")) {
        cout << cmdline_options << std::endl;
        return false;
    }

    if (vm.count("help")) {
        cout << visible << std::endl;
        return false;
    }

    if (vm.count("version")) {
        cout << bformat("S-W Graph Constructor, part of IgReC version %s; git version: %s") % build_info::version % build_info::git_hash7 << std::endl;
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
            store(po::command_line_parser(argc, argv).
                  options(cmdline_options).positional(p).run(), vm);
        }
    }

    try {
        notify(vm);
    } catch (po::error &e) {
        cout << "Parser error: " << e.what() << std::endl;
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

    SWGCParam args;
    try {
        if (!parse_cmd_line_arguments(argc, argv, args)) {
            return 0;
        }
    } catch(std::exception& e) {
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

    INFO("Read length checking");
    size_t required_read_length = (args.strategy != 0) ? (args.k * (args.tau + args.strategy)) : 0;
    size_t required_read_length_for_single_strategy = args.k * (args.tau + 1);
    size_t required_read_length_for_double_strategy = args.k * (args.tau + 2);

    size_t discarded_reads = 0;
    size_t discarded_reads_single = 0;
    size_t discarded_reads_double = 0;
    for (const auto &read : input_reads) {
        discarded_reads += length(read) < required_read_length;
        discarded_reads_single += length(read) < required_read_length_for_single_strategy;
        discarded_reads_double += length(read) < required_read_length_for_double_strategy;
    }

    int saved_reads_single = static_cast<int>(discarded_reads) - static_cast<int>(discarded_reads_single);
    int saved_reads_double = static_cast<int>(discarded_reads) - static_cast<int>(discarded_reads_double);

    if (saved_reads_single > 0.05 * static_cast<double>(input_reads.size())) {
        if (saved_reads_single - saved_reads_double < 0.05 * static_cast<double>(input_reads.size())) {
            INFO(bformat("Choosing <<double>> strategy for saving %d reads")
                 % saved_reads_double);
            args.strategy = 2;
            discarded_reads = discarded_reads_double;
        } else {
            INFO(bformat("Choosing <<single>> strategy for saving %d reads")
                 % saved_reads_single);
            args.strategy = 1;
            discarded_reads = discarded_reads_single;
        }
    }

    if (discarded_reads) {
        WARN(bformat("Discarded reads %d") % discarded_reads);
    }

    omp_set_num_threads(args.nthreads);
    INFO(bformat("Truncated distance graph construction using %d threads starts") % args.nthreads);
    INFO("Construction of candidates graph");

    INFO("Strategy " << args.strategy << " was chosen");

    auto dist_fun = [&args](const Dna5String& s1, const Dna5String& s2) -> unsigned {
        auto delta = [&args](int l) -> int { return (bool)(l)*2 * args.tau; };
        auto lizard_tail = [&args, &delta](int l) -> int { return args.ignore_tails ? 0 : -delta(l); };
        return -half_sw_banded(s1, s2, 0, -1, -1, lizard_tail, args.max_indels);
    };

    if (args.reference_file == "") {
        INFO("K-mer index construction");
        auto kmer2reads = kmerIndexConstruction(input_reads, args.k);

        size_t num_of_dist_computations;
        auto dist_graph = tauDistGraph(input_reads,
                                       kmer2reads,
                                       dist_fun,
                                       args.tau, args.k,
                                       args.strategy,
                                       num_of_dist_computations);

        INFO("Simularity computations: " << num_of_dist_computations << ", average " << \
             static_cast<double>(num_of_dist_computations) / static_cast<double>(input_reads.size()) << " per read");

        size_t num_of_edges = numEdges(dist_graph);
        INFO("Edges found: " << num_of_edges);
        INFO("Strategy efficiency: " << static_cast<double> (num_of_edges) / static_cast<double>(num_of_dist_computations));

        // Output
        if (args.export_abundances) {
            INFO("Saving graph (with abundances)");
            auto abundances = find_abundances(input_ids);
            write_metis_graph(dist_graph, abundances, args.output_file);
        } else {
            INFO("Saving graph (without abundances)");
            write_metis_graph(dist_graph, args.output_file);
        }
    } else {
        SeqFileIn seqFileIn_reference(args.reference_file.c_str());
        std::vector<CharString> reference_ids;
        std::vector<Dna5String> reference_reads;

        INFO("Reading input reads starts");
        readRecords(reference_ids, reference_reads, seqFileIn_reference);
        INFO(reference_reads.size() << " reads were extracted from " << args.reference_file);

        INFO("K-mer index construction");
        auto kmer2reads = kmerIndexConstruction(reference_reads, args.k);

        size_t num_of_dist_computations;
        auto dist_graph = tauMatchGraph(input_reads,
                                        reference_reads,
                                        kmer2reads,
                                        dist_fun,
                                        args.tau, args.k,
                                        args.strategy,
                                        num_of_dist_computations);

        INFO("Simularity computations: " << num_of_dist_computations << ", average " << \
             static_cast<double>(num_of_dist_computations) / static_cast<double>(input_reads.size()) << " per read");

        size_t num_of_edges = numEdges(dist_graph, false);
        INFO("Edges found: " << num_of_edges);
        INFO("Strategy efficiency: " << static_cast<double> (num_of_edges) / static_cast<double>(num_of_dist_computations));

        // Output
        if (args.export_abundances) {
            INFO("Saving graph (with abundances)");
            auto abundances = find_abundances(input_ids);
            write_metis_graph(dist_graph, abundances, args.output_file, false);
        } else {
            INFO("Saving graph (without abundances)");
            write_metis_graph(dist_graph, args.output_file, false);
        }
    }

    INFO("Graph was written to " << args.output_file);

    INFO("Running time: " << running_time_format(pc));

    return 0;
}

// vim: ts=4:sw=4
