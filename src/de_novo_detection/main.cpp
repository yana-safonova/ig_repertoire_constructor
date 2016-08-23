
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <logger/logger.hpp>
#include <logger/log_writers.hpp>

// #include "construction.hpp"
#include "partial_order_graph.hpp"

namespace po = boost::program_options;

struct options {
    std::string input_file;
    std::string graph_prefix;
    std::string output_directory;
    std::string logging_level;
    bool no_plots;
    bool show_sequences;

    bool use_forward_reads;
    bool use_rc_reads;

    size_t k_step;
    size_t k_stop;
};


void parse_command_line(int argc, char** argv, options& opts) {
    std::string name = path::basename(argv[0]);
    std::string usage =
    R"(De novo IGHV construction

Usage:
  )" + name + R"( -i <in.*> [-g <graph>] -o <out-dir> [options]
)";

    po::options_description required_options("Required options");
    required_options.add_options()
        (",i", po::value<std::string>(&opts.input_file)->value_name("<in.*>"), "Input sequences")
        (",g", po::value<std::string>(&opts.graph_prefix)->value_name("<graph>"),
               "Graph prefix. Should exist files \"<graph>.csv\" and \"<graph>.dot\"")
        (",o", po::value<std::string>(&opts.output_directory)->value_name("<out-dir>"), "Output directory")
        ;

    pog::pog_parameters& parameters = pog::pog_parameters::instance();
    parameters.mismatch_penalty = -1.4f;
    parameters.gap_penalty = -2.0f;
    std::string orientation = "f";
    parameters.bulges_hamming_ratio = .3f;
    parameters.coverage_threshold = 10.0f;

    po::options_description pog_options("Partial order graph options");
    pog_options.add_options()
        ("mismatch-penalty", po::value<float>(&parameters.mismatch_penalty)->value_name("FLOAT"),
            "Mismatch penalty   (default: -1.4)")
        ("gap-penalty", po::value<float>(&parameters.gap_penalty)->value_name("FLOAT"),
            "Gap penalty        (default:   -2)")
        ("orientation", po::value<std::string>(&orientation)->value_name("{f,r,fr}"),
            "Reads orientation, {f,r,fr}: forward, reverse complement or both, (default: fr")
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
        ("log-level", po::value<std::string>(&opts.logging_level)->value_name("LEVEL")->default_value("info"),
            "Logging level {trace, debug, info, warn, error} (default: info)")
        ;

    po::options_description full_description;
    full_description
        .add(required_options)
        .add(pog_options)
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
    if (!path::check_existence(opts.output_directory) && !path::make_dir(opts.output_directory)) {
        std::cerr << "Could not create directory " << opts.output_directory << std::endl;
        exit(1);
    }
    opts.output_directory += '/';

    opts.no_plots = vm.count("no-plots");
    opts.show_sequences = vm.count("show-sequences");

    opts.use_forward_reads = orientation.find('f') != std::string::npos;
    opts.use_rc_reads = orientation.find('r') != std::string::npos;
}

void construct(options& opts) {
    auto& params = pog::pog_parameters::instance();

    INFO("Constructing de novo")
    INFO("k = 1")
    params.set_kmer_size(1);

    pog::partial_order_graph graph;
    if (opts.graph_prefix.empty())
        pog::from_file(graph, opts.input_file, opts.use_forward_reads, opts.use_rc_reads);
    else
        graph.load(opts.graph_prefix);

    graph.remove_low_covered();
    pog::save_graph(graph, opts.output_directory + "1", opts.show_sequences);
    if (!opts.no_plots)
        pog::draw_graph(opts.output_directory + "1");
    pog::correct_reads(graph, opts.input_file, opts.output_directory + "1.fa",
                       opts.use_forward_reads,
                       opts.use_rc_reads);
    // graph.clear();
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

int main(int argc, char** argv) {
    options opts;
    parse_command_line(argc, argv, opts);
    create_console_logger(string_to_level(opts.logging_level));
    construct(opts);
}
