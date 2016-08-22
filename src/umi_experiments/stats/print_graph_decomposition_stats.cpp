#include "../../graph_utils/graph_io.hpp"
#include "graph_stats.hpp"
#include "../../ig_tools/utils/string_tools.hpp"
#include "utils.hpp"
#include <segfault_handler.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

bool readArgs(int argc, char **argv, std::string& reads_file, std::string& graph_file, std::string& output_file, size_t& size_to_print) {
    namespace po = boost::program_options;
    po::options_description cmdl_options("Is this needed?");
    cmdl_options.add_options()
            ("help,h", "print help message")
            ("reads,r", po::value<std::string>(&reads_file)->required(), "input file with reads")
            ("graph,g", po::value<std::string>(&graph_file)->required(), "input file with graph")
            ("output,o", po::value<std::string>(&output_file)->default_value(""), "output file with stats")
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

void extract_abundances(const std::vector<seqan::CharString>& ids, std::vector<size_t>& abundances) {
    // error somewhere here on age3
    for (auto& id : ids) {
        const std::string& id_string = seqan_string_to_string(id);
        const std::string& abundance_string = id_string.substr(id_string.find_last_of("___") + 1);
        abundances.push_back(stoull(abundance_string));
    }
}

int main(int argc, char **argv) {
    segfault_handler sh;
    create_console_logger();
    std::string reads_file;
    std::string graph_file;
    std::string output_file;
    size_t size_to_print;
    try {
        if (!readArgs(argc, argv, reads_file, graph_file, output_file, size_to_print)) {
            return 0;
        }
    } catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    INFO("Reading fastq");
    seqan::SeqFileIn seqFileIn_input(reads_file.c_str());
    std::vector<seqan::CharString> input_ids;
    std::vector<seqan::Dna5String> input_reads;
    readRecords(input_ids, input_reads, seqFileIn_input);
    std::vector<size_t> abundances;
    extract_abundances(input_ids, abundances);

    INFO("Reading graph")
    const SparseGraphPtr graph = GraphReader(graph_file).CreateGraph();

    INFO("Analyzing structure");
    const GraphStats stats = GraphStats::GetStats(input_reads, abundances, graph);

    if (output_file != "") {
        std::ofstream output(output_file);
        output << stats.ToString(size_to_print);
        INFO("Stats printed to " << output_file);
    } else {
        INFO(stats.ToString(size_to_print));
    }
}
