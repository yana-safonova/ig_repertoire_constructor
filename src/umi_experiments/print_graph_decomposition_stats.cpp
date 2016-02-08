#include <logger/log_writers.hpp>
#include "../graph_utils/graph_io.hpp"
#include "graph_stats.hpp"

bool readArgs(int argc, char **argv, string& input_file, size_t& size_to_print) {
    if (argc < 2 || argc > 3) {
        cout << "Usage: <input file> [<size of unclassified components to print>]" << endl;
        return false;
    }
    input_file = argv[1];
    size_to_print = argc == 3 ? atoi(argv[2]) : 0;
    return true;
}

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char **argv) {
    create_console_logger();
    string input_file = "";
    size_t size_to_print = 0;
    if (!readArgs(argc, argv, input_file, size_to_print)) {
        return 0;
    }
    const SparseGraphPtr graph = GraphReader(input_file).CreateGraph();
    INFO("Searching and removing stars");
    const GraphStats stats = GraphStats::GetStats(graph);

    INFO(stats.ToString(size_to_print));
}
