#include <logger/log_writers.hpp>
#include "../graph_utils/graph_io.hpp"
#include "graph_stats.hpp"
#include <segfault_handler.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

bool readArgs(int argc, char **argv, string& input_file, string& output_file, size_t& size_to_print) {
    po::options_description cmdl_options("Is this needed?");
    cmdl_options.add_options()
            ("help,h", "print help message")
            ("input,i", po::value<string>(&input_file)->required(), "input file with graph")
            ("output,o", po::value<string>(&output_file)->default_value(""), "output file with stats")
            ("size,s", po::value<size_t>(&size_to_print)->default_value(std::numeric_limits<size_t>::max()), "minimal size of printed component")
            ;
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(cmdl_options).run(), vm);
    if (vm.count("help") || argc == 1) {
        cout << cmdl_options << endl;
        return false;
    }
    po::notify(vm);
    return true;
}

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char **argv) {
    segfault_handler sh;
    create_console_logger();
    string input_file;
    string output_file;
    size_t size_to_print;
    try {
        if (!readArgs(argc, argv, input_file, output_file, size_to_print)) {
            return 0;
        }
    } catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    const SparseGraphPtr graph = GraphReader(input_file).CreateGraph();
    INFO("Searching and removing stars");
    const GraphStats stats = GraphStats::GetStats(graph);

    if (output_file != "") {
        std::ofstream output(output_file);
        output << stats.ToString(size_to_print);
        INFO("Stats printed to " << output_file);
    } else {
        INFO(stats.ToString(size_to_print));
    }
}
