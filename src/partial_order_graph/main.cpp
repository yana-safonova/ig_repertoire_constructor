#include <iostream>
#include <string>
#include <vector>
#include <seqan/seq_io.h>
#include <boost/program_options.hpp>
#include <logger/logger.hpp>
#include <logger/log_writers.hpp>
#include <copy_file.hpp>
#include <stdlib.h>

#include "pog_parameters.hpp"
#include "partial_order_graph.hpp"

namespace po = boost::program_options;

struct options {
    std::string input_file;
    std::string directory;
    std::string logging_level;
    std::string prefix;

    bool no_plots;
    bool show_sequences;
    bool shrink_upaths;
    bool shrink_bulges;
    bool correct;
    bool remove_low_covered;

    bool use_forward_reads;
    bool use_rc_reads;
};

void parse_command_line(int argc, char** argv, options& opts) {
    std::string name = path::basename(argv[0]);
    std::string usage =
    R"(Partial order graph construction and shrinkage

Usage:
  )" + name + R"( -i <in.*> [-d <directory> | -p <prefix>] -k INT [options]
)";

    size_t kmer_size = -1;
    po::options_description required_options("Required options");
    required_options.add_options()
        (",i", po::value<std::string>(&opts.input_file)->value_name("<in.*>"), "Input sequences")
        (",d", po::value<std::string>(&opts.directory)->value_name("<directory>"), "Output directory")
        (",p", po::value<std::string>(&opts.prefix)->value_name("<prefix>"), "Prefix for output files")
        (",k", po::value<size_t>(&kmer_size)->value_name("INT"), "k-mer size (default: 10)")
        ;

    pog::pog_parameters& parameters = pog::pog_parameters::instance();

    parameters.mismatch_penalty = -1.4f;
    parameters.gap_penalty = -2.0f;
    std::string orientation = "f";
    po::options_description construction_options("Construction options");
    construction_options.add_options()
        ("mismatch-penalty", po::value<float>(&parameters.mismatch_penalty)->value_name("FLOAT"),
            "Mismatch penalty   (default: -1.4)")
        ("gap-penalty", po::value<float>(&parameters.gap_penalty)->value_name("FLOAT"),
            "Gap penalty        (default:   -2)")
        ("orientation", po::value<std::string>(&orientation)->value_name("{f,r,fr}"),
            "Reads orientation, {f,r,fr}: forward, reverse complement or both, (default: f)")
        ;

    parameters.bulges_hamming_ratio = .3f;
    parameters.coverage_threshold = 10.0f;
    po::options_description shrinkage_options("Shrinkage options");
    shrinkage_options.add_options()
        ("shrink-upaths", "Shrink unambigous paths")
        ("shrink-bulges", "Shrink bulges")
        ("remove-low-covered", "Remove low covered nodes")
        ("correct", "Correct reads")
        ("bulges-hamming", po::value<float>(&parameters.bulges_hamming_ratio)->value_name("FLOAT"),
            "Hamming distance between nodes <= ceil(bulges:hamming * length): condition "
            "to join nodes during bulge shrinking (default: 0.3)")
        ("cov-threshold", po::value<float>(&parameters.coverage_threshold)->value_name("FLOAT"),
            "Coverage threshold to remove vertex during bulge shrinking (default: 10)")
        ;

    po::options_description output_options("Output options");
    output_options.add_options()
        ("no-plots", "Do not call graphviz")
        ("show-sequences", "Show sequences on the graph plot")
        ;

    po::options_description other_options("Other");
    other_options.add_options()
        ("help,h", "Show this message")
        ("log-level", po::value<std::string>(&opts.logging_level)->value_name("LEVEL"),
            "Logging level {trace, debug, info, warn, error} (default: info)")
        ;

    po::options_description full_description;
    full_description
        .add(required_options)
        .add(construction_options)
        .add(shrinkage_options)
        .add(output_options)
        .add(other_options)
        ;

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, full_description), vm);
        po::notify(vm);
    } catch (const po::error &e) {
        std::cerr << "Command-line parser error: " << e.what() << std::endl;
        exit(1);
    } catch (const std::exception &e) {
        std::cerr << "Unknown exception: " << e.what() << std::endl;
        exit(1);
    }

    if (vm.count("help")) {
        std::cout << usage << full_description << "\n";
        exit(0);
    }

    if (kmer_size == static_cast<size_t>(-1)) {
        std::cerr << "K-mer size \"-k INT\" is required" << std::endl;
        exit(1);
    }
    parameters.set_kmer_size(kmer_size);

    if (opts.input_file.empty()) {
        std::cerr << "Input sequences \"-i <in.*>\" are required" << std::endl;
        exit(1);
    }
    if (opts.directory.empty() == opts.prefix.empty()) {
        std::cerr << "Specify either \"-d <dir>\" or \"-p <prefix>\"" << std::endl;
        exit(1);
    }

    if (!opts.directory.empty()) {
        if (!path::check_existence(opts.directory) && !path::make_dir(opts.directory)) {
            std::cerr << "Could not create directory " << opts.directory << std::endl;
            exit(1);
        }
        opts.prefix = opts.directory + "/" + std::to_string(kmer_size);
    }

    opts.no_plots = vm.count("no-plots");
    opts.show_sequences = vm.count("show-sequences");
    opts.shrink_upaths = vm.count("shrink-upaths");
    opts.shrink_bulges = vm.count("shrink-bulges");
    opts.remove_low_covered = vm.count("remove-low-covered");
    opts.correct = vm.count("correct");

    opts.use_forward_reads = orientation.find('f') != std::string::npos;
    opts.use_rc_reads = orientation.find('r') != std::string::npos;
}

logging::level string_to_level(std::string const& level_str) {
    using namespace logging;
    if (level_str == "" || level_str == "info")
        return L_INFO;
    if (level_str == "trace")
        return L_TRACE;
    if (level_str == "debug")
        return L_DEBUG;
    if (level_str == "warn")
        return L_WARN;
    if (level_str == "error")
        return L_ERROR;
    std::cerr << "Unsupported logging level: " << level_str << std::endl;
    exit(1);
}

void create_console_logger(logging::level log_level) {
    logging::logger *lg = logging::create_logger("", log_level);
    lg->add_writer(std::make_shared<logging::console_writer>());
    logging::attach_logger(lg);
}

void correct_reads(options& opts) {
    pog::partial_order_graph graph;
    pog::from_file(graph, opts.input_file, opts.use_forward_reads, opts.use_rc_reads);
    INFO("Graph created. Nodes: " << graph.nodes_count());

    graph.remove_low_covered();
    pog::save_graph(graph, opts.prefix, opts.show_sequences);
    if (!opts.no_plots)
        pog::draw_graph(opts.prefix);

    pog::correct_reads(graph, opts.input_file, opts.prefix + "corr.fa",
                       opts.use_forward_reads,
                       opts.use_rc_reads);
}

int main(int argc, char** argv) {
    options opts;
    parse_command_line(argc, argv, opts);

    create_console_logger(string_to_level(opts.logging_level));
    if (opts.correct) {
        correct_reads(opts);
        return 0;
    }

    pog::partial_order_graph graph;
    pog::from_file(graph, opts.input_file, opts.use_forward_reads, opts.use_rc_reads);
    INFO("Graph created. Nodes: " << graph.nodes_count());


    if (opts.remove_low_covered)
        graph.remove_low_covered();
    if (opts.shrink_upaths)
        graph.shrink_upaths();
    if (opts.shrink_bulges)
        graph.shrink_bulges();

    pog::save_graph(graph, opts.prefix, opts.show_sequences);
    if (!opts.no_plots)
        pog::draw_graph(opts.prefix);

    return 0;
}
