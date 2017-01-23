//
// Created by Andrew Bzikadze on 5/27/16.
//

#include <boost/program_options.hpp>
#include <logger/logger.hpp>
#include "command_line_routines.hpp"
#include "shm_config.hpp"

namespace shm_kmer_matrix_estimator {

bool command_line_requires_parsing(int argc, char **argv) {
    if (argc == 1)
        return false;
    if (argc > 2)
        return true;
    return std::string(argv[1]) == "--help" or std::string(argv[1]) == "-h";
}

// cfg contains default values from config file
void parse_command_line_args(shm_config &cfg, int argc, char **argv) {
    if (!command_line_requires_parsing(argc, argv))
        return;

    namespace po = boost::program_options;

    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
        ("help,h", "produce help message")
        ("v_alignments,v", po::value<std::string>(&cfg.io.input.v_alignments)->
             default_value(cfg.io.input.v_alignments),
         "name of an input v_alignments file (FASTA)")
        ("cdr_details,d", po::value<std::string>(&cfg.io.input.cdr_details)->
             default_value(cfg.io.input.cdr_details),
         "name of an input cdr details file (FASTA)")
        ("output-filename-fr,f", po::value<std::string>(&cfg.io.output.output_filename_fr)->
             default_value(cfg.io.output.output_filename_fr),
         "output filename for FR statistics")
        ("output-filename-cdr,c", po::value<std::string>(&cfg.io.output.output_filename_cdr)->
             default_value(cfg.io.output.output_filename_cdr),
         "output filename for CDR statistics")
        ("mutation-strategy,s", po::value<shm_config::mutations_strategy_params::MutationsStrategyMethod>
             (&cfg.mfp.mutations_strategy_method)->
             default_value(cfg.mfp.mutations_strategy_method,
                           cfg.mfp.mutation_strategy_method_names[
                               static_cast<size_t>(cfg.mfp.mutations_strategy_method)]),
         "mutation strategy: Trivial or NoKNeighbours");
    po::options_description cmdline_options("All command line options");
    cmdline_options.add(generic);

    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);

    if (vm.count("help")) {
        std::cout << generic;
        exit(0);
    }

    try {
        po::notify(vm);
    } catch (const po::error &e) {
        ERROR("Command-line parser error: " << e.what());
        exit(1);
    } catch (const std::exception &e) {
        ERROR("Unknown exception: " << e.what());
        exit(1);
    }
}

} // End namespace shm_kmer_matrix_estimator
