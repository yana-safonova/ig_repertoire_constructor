#include <perfcounter.hpp>
#include <segfault_handler.hpp>
#include <logger/log_writers.hpp>
#include <boost/program_options.hpp>
#include "error_analyzer.hpp"

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

void parse_options(int argc, const char* const* argv, std::string& input_file_path, std::string& error_stats_dir_path) {
    namespace po = boost::program_options;
    po::options_description cmdline_options("Command line options");
    cmdline_options.add_options()
            ("input-file,i", po::value<std::string>(&input_file_path)->required(),
             "name of the input file with reads (should be cleaned by vj-finder)")
            ("output-dir,o", po::value<std::string>(&error_stats_dir_path),
             "if present, error stats will be output to this directory");
    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
    po::notify(vm);
}

int main(int argc, const char* const* argv) {
    perf_counter pc;
    segfault_handler sh;
    create_console_logger();

    std::string input_file_path;
    std::string error_stats_dir_path;
    parse_options(argc, argv, input_file_path, error_stats_dir_path);

    ErrorAnalyzer error_analyzer;
    error_analyzer.readData(input_file_path);
    error_analyzer.performAnalysis(error_stats_dir_path);
}