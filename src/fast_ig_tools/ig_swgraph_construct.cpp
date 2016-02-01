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


template<typename T, typename Tf>
Graph tauDistGraph(const std::vector<T> &input_reads,
                   const KmerIndex &kmer2reads,
                   const Tf &dist_fun,
                   int tau,
                   int K,
                   Strategy strategy,
                   size_t &num_of_dist_computations) {
    Graph g(input_reads.size());

    std::atomic<size_t> atomic_num_of_dist_computations;
    atomic_num_of_dist_computations = 0;

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 8))
    for (size_t j = 0; j < input_reads.size(); ++j) {
        auto cand = find_candidates(input_reads[j], kmer2reads, input_reads.size(), tau, K, strategy);

        size_t len_j = length(input_reads[j]);

        for (size_t i : cand) {
            size_t len_i = length(input_reads[i]);
            if (len_j < len_i || (len_i == len_j && j < i)) {
                int dist = dist_fun(input_reads[j], input_reads[i]);

                atomic_num_of_dist_computations += 1;

                if (dist <= tau) {
                    g[j].push_back( { i, dist } );
                }
            }
        }
    }

    // Undirecting
    auto gg = g;
    for (size_t i = 0; i < gg.size(); ++i) {
        for (const auto &_ : gg[i]) {
            g[_.first].push_back( { i, _.second } );
        }
    }
    gg.clear(); // Free memory

    SEQAN_OMP_PRAGMA(parallel for schedule(guided, 8))
    for (size_t j = 0; j < g.size(); ++j) {
        remove_duplicates(g[j]);
    }

    num_of_dist_computations = atomic_num_of_dist_computations;

    return g;
}

bool parse_cmd_line_arguments(int argc, char **argv, std::string &input_file, std::string &output_file, int &K, int &tau,
                              int &nthreads, int &strategy_int, int &max_indels, bool &export_abundances) {
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
             "file for outputted truncated dist-graph in METIS format")
            ("export-abundances,A", "export read abundances to output graph file")
            ("no-export-abundances", "don't export read abundances to output graph file (default)")
            ;

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
            ("word-size,k", po::value<int>(&K)->default_value(K),
             "word size for k-mer index construction")
            ("strategy,S", po::value<int>(&strategy_int)->default_value(strategy_int),
             "strategy type (0 --- naive, 1 --- single, 2 --- pair, 3 --- triple)")
            ("tau", po::value<int>(&tau)->default_value(tau),
             "maximum distance value for truncated dist-graph construction")
            ("max-indels", po::value<int>(&max_indels)->default_value(max_indels),
             "maximum number of indels in Levenshtein distance")
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

    if (vm.count("export-abundances")) {
        export_abundances = true;
    }

    if (vm.count("no-export-abundances")) {
        export_abundances = false;
    }
    return true;
}

int main(int argc, char **argv) {
    segfault_handler sh;
    perf_counter pc;
    create_console_logger("");

    INFO("Command line: " << join_cmd_line(argc, argv));

    int K = 10; // anchor length
    int tau = 4;
    int nthreads = 4;
    std::string input_file = "cropped.fa";
    std::string output_file = "output.graph";
    int strategy_int = Strategy::TRIPLE;
    int max_indels = 0;
    bool export_abundances = false;

    try {
        if (!parse_cmd_line_arguments(argc, argv, input_file, output_file, K, tau, nthreads, strategy_int, max_indels, export_abundances)) {
            return 0;
        }
    } catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    INFO("Input reads: " << input_file);

    INFO("K = " << K << ", tau = " << tau);

    SeqFileIn seqFileIn_input(input_file.c_str());
    std::vector<CharString> input_ids;
    std::vector<Dna5String> input_reads;

    INFO("Reading input reads starts");
    readRecords(input_ids, input_reads, seqFileIn_input);
    INFO(input_reads.size() << " reads were extracted from " << input_file);

    INFO("Read length checking");
    size_t required_read_length = (strategy_int != 0) ? (K * (tau + strategy_int)) : 0;
    size_t required_read_length_for_single_strategy = K * (tau + 1);
    size_t required_read_length_for_double_strategy = K * (tau + 2);

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
            strategy_int = 2;
            discarded_reads = discarded_reads_double;
        } else {
            INFO(bformat("Choosing <<single>> strategy for saving %d reads")
                 % saved_reads_single);
            strategy_int = 1;
            discarded_reads = discarded_reads_single;
        }
    }

    if (discarded_reads) {
        WARN(bformat("Discarded reads %d") % discarded_reads);
    }

    INFO("K-mer index construction");
    auto kmer2reads = kmerIndexConstruction(input_reads, K);

    omp_set_num_threads(nthreads);
    INFO(bformat("Truncated distance graph construction using %d threads starts") % nthreads);
    INFO("Construction of candidates graph");

    Strategy strategy = Strategy(strategy_int);
    INFO(toCString(strategy) << " was chosen");

    auto dist_fun = [max_indels](const Dna5String &s1, const Dna5String &s2) -> int {
        auto lizard_tail = [](int l) -> int { return 0*l; };
        return -half_sw_banded(s1, s2, 0, -1, -1, lizard_tail, max_indels);
    };

    size_t num_of_dist_computations;
    auto dist_graph = tauDistGraph(input_reads,
                                   kmer2reads,
                                   dist_fun,
                                   tau, K,
                                   strategy,
                                   num_of_dist_computations);


    INFO("Simularity computations: " << num_of_dist_computations << ", average " << \
         static_cast<double>(num_of_dist_computations) / input_reads.size() << " per read");

    size_t num_of_edges = numEdges(dist_graph);
    INFO("Edges found: " << num_of_edges);
    INFO("Strategy efficiency: " << static_cast<double> (num_of_edges) / num_of_dist_computations);

    // Output
    if (export_abundances) {
        INFO("Saving graph (with abundances)");
        auto abundances = find_abundances(input_ids);
        write_metis_graph(dist_graph, abundances, output_file);
    } else {
        INFO("Saving graph (without abundances)");
        write_metis_graph(dist_graph, output_file);
    }
    INFO("Graph was written to " << output_file);

    INFO("Running time: " << running_time_format(pc));

    return 0;
}

// vim: ts=4:sw=4
