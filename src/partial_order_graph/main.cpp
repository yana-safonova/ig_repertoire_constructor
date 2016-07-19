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
    std::string output_directory;
    std::string logging_level;
    std::string prefix;
    bool no_plots;
    bool show_sequences;
    bool wo_reverse_complement;
    bool shrink_upaths;
    bool shrink_bulges;
    bool remove_low_covered;

    bool use_forward_reads;
    bool use_rc_reads;
};

void parse_command_line(int argc, char** argv, options& opts) {
    std::string name = path::basename(argv[0]);
    std::string usage =
    R"(Partial order graph construction and shrinkage

Usage:
  )" + name + R"( -i <in.*> -o <out-dir> [options]
)";

    po::options_description required_options("Required options");
    required_options.add_options()
        (",i", po::value<std::string>(&opts.input_file)->value_name("<in.*>"), "Input sequences")
        (",o", po::value<std::string>(&opts.output_directory)->value_name("<out-dir>"), "Output directory")
        ;

    pog::pog_parameters& parameters = pog::pog_parameters::instance();
    size_t kmer_size = parameters.get_kmer_size();

    std::string orientation = "fr";
    po::options_description construction_options("Construction options");
    construction_options.add_options()
        (",k", po::value<size_t>(&kmer_size)->value_name("INT"),
            "k-mer size         (default:   10)")
        ("mismatch-penalty", po::value<float>(&parameters.mismatch_penalty)->value_name("FLOAT"),
            "Mismatch penalty   (default: -1.4)")
        ("gap-penalty", po::value<float>(&parameters.gap_penalty)->value_name("FLOAT"),
            "Gap penalty        (default:   -2)")
        ("orientation", po::value<std::string>(&orientation)->value_name("{f,r,fr}"),
            "Reads orientation, {f,r,fr}: forward, reverse complement or both, (default: fr)")
        ;

    po::options_description shrinkage_options("Shrinkage options");
    shrinkage_options.add_options()
        ("shrink-upaths", "Shrink unambigous paths")
        ("shrink-bulges", "Shrink bulges")
        ("remove-low-covered", "Remove low covered nodes")
        ("bulges-hamming", po::value<float>(&parameters.bulges_hamming_ratio)->value_name("FLOAT"),
            "Hamming distance between nodes <= ceil(bulges:hamming * length): condition "
            "to join nodes during bulge shrinking (default: 0.3)")
        ("coverage-threshold", po::value<float>(&parameters.coverage_threshold)->value_name("FLOAT"),
            "Coverage threshold to remove vertex during bulge shrinking (default: 10)")
        ;

    po::options_description output_options("Output options");
    output_options.add_options()
        ("no-plots", "Do not call graphviz")
        ("show-sequences", "Show sequences on the graph plot")
        ("prefix,p", po::value<std::string>(&opts.prefix)->value_name("PREFIX"), "File prefix")
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

    if (opts.input_file.empty()) {
        std::cerr << "Input sequences \"-i <in.*>\" are required" << std::endl;
        exit(1);
    }
    if (opts.output_directory.empty()) {
        std::cerr << "Output directory \"-o <out-dir>\" is required" << std::endl;
        exit(1);
    }

    parameters.set_kmer_size(kmer_size);
    opts.no_plots = vm.count("no-plots");
    opts.show_sequences = vm.count("show-sequences");
    opts.wo_reverse_complement = vm.count("wo-rc");
    opts.shrink_upaths = vm.count("shrink-upaths");
    opts.shrink_bulges = vm.count("shrink-bulges");
    opts.remove_low_covered = vm.count("remove-low-covered");

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

void plot_if_necessary(options& opts) {
    if (opts.no_plots)
        return;
    INFO("Drawing plot " + opts.prefix + ".pdf");
    system(("dot " + opts.output_directory + opts.prefix + ".dot -T pdf -o "
                + opts.output_directory + opts.prefix + ".pdf").c_str());
}

int main(int argc, char** argv) {
    options opts;
    parse_command_line(argc, argv, opts);
    opts.output_directory += "/";
    if (opts.prefix != "")
        opts.prefix += ".";
    opts.prefix += std::to_string(pog::pog_parameters::instance().get_kmer_size());

    if (!path::check_existence(opts.output_directory) && !path::make_dir(opts.output_directory)) {
        std::cerr << "Could not create directory " << opts.output_directory << std::endl;
        return 1;
    }

    create_console_logger(string_to_level(opts.logging_level));
    pog::partial_order_graph graph = pog::from_file(opts.input_file, opts.use_forward_reads, opts.use_rc_reads);
    INFO("Graph created. Nodes: " << graph.nodes_count());

    if (opts.shrink_upaths)
        graph.shrink_upaths();
    if (opts.shrink_bulges)
        graph.shrink_bulges();
    if (opts.remove_low_covered)
        graph.remove_low_covered();

    INFO("Saving");
    graph.save_nodes(opts.output_directory + opts.prefix + ".csv");
    graph.save_dot(opts.output_directory + opts.prefix, opts.show_sequences);
    plot_if_necessary(opts);
    return 0;
}
