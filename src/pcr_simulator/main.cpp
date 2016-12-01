#include <logger/logger.hpp>
#include <logger/log_writers.hpp>
#include <segfault_handler.hpp>
#include <perfcounter.hpp>
#include <boost/program_options.hpp>
#include "options.hpp"
#include "pcr_simulator.hpp"

Options parse_options(int argc, const char* const* argv) {
    Options options;
    Options::SimulationOptions& simulationOptions = options.simulation_options;
    namespace po = boost::program_options;
    po::options_description cmdline_options("Command line options");
    cmdline_options.add_options()
            ("input-file,i", po::value<std::string>(&options.repertoire_file_path)->required(),
             "name of the input repertoire file (FASTA|FASTQ)")
            ("output-dir,o", po::value<std::string>(&options.output_dir_path)->required(),
             "path to the output directory")
            ("umi-length,l", po::value<size_t>(&simulationOptions.barcode_length)->default_value(15),
             "length of generated barcodes (defaults to 14)")
            ("pcr-cycles,c", po::value<size_t>(&simulationOptions.cycles_count)->default_value(25),
             "number of PCR cycles to simulate (defaults to 25)")
            ("pcr-error1,e", po::value<double>(&simulationOptions.error_prob_first)->required(),
             "probability of PCR error on the first PCR cycle")
            ("pcr-error2,E", po::value<double>(&simulationOptions.error_prob_last)->required(),
             "probability of PCR error on the last PCR cycle (probability on the other cycles are interpolated)")
            ("pcr-rate,r", po::value<double>(&simulationOptions.amplification_rate)->default_value(0.1),
             "probability for each molecule to be amplified on each PCR cycle (defaults to 0.1)")
            ("chimeras-rate,k", po::value<double>(&simulationOptions.chimeras_rate)->default_value(0.001),
             "share of chimeric reads appering on each cycle (defaults to 0.001)")
            ("output-limit,m", po::value<size_t>(&options.output_estimation_limit)->default_value(static_cast<size_t>(1e8)),
             "the program will exit if expected number of reads exceeds this parameter (defaults to 100,000,000)")
            ("barcode-position,b", po::value<size_t>(&simulationOptions.barcode_position)->default_value(3),
             "indicator of barcode position in the read, used for chimeras simulation. "
             "1 for barcode going with the left half, 2 for the right half, 3 for random choice (defaults to 3)")
            ;
    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
    po::notify(vm);
    return options;
}

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, const char* const* argv) {
    perf_counter pc;
    segfault_handler sh;
    create_console_logger();
    std::srand(495634);

    const Options& options = parse_options(argc, argv);

    PcrSimulator simulator(options.simulation_options);
    simulator.ReadRepertoire(options.repertoire_file_path);
    simulator.Amplify(options.output_estimation_limit);
    simulator.WriteResults(options.output_dir_path);
}
