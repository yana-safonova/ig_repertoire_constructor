#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#include "fast_ig_tools.hpp"

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::CharString;



struct IgQuastParameters {
    int K = 36; // anchor length
    int tau = 3;
    int nthreads = 4;
    std::string input_file1 = "input1.fa";
    std::string input_file2 = "input2.fa";
    std::string output_file1 = "output1.match";
    std::string output_file2 = "output2.match";
    int strategy_int = Strategy::TRIPLE;

int parse_cmd_line 


    // Parse cmd-line arguments
    try {
        std::string config_file = "";

        // Declare a group of options that will be
        // allowed only on command line
        po::options_description generic("Generic options");
        generic.add_options()
            ("version,v", "print version string")
            ("help,h", "produce help message")
            ("config,c", po::value<std::string>(&config_file)->default_value(config_file),
             "name of a file of a configuration")
            ("input-file1,i", po::value<std::string>(&input_file1),
             "name of the first input file (FASTA|FASTQ)")
            ("input-file2,I", po::value<std::string>(&input_file2),
             "name of the second input file (FASTA|FASTQ)")
            ("output-file1,o", po::value<std::string>(&output_file1)->default_value(output_file1),
             "file for outputted match for reads of the first file to reads of the second one")
            ("output-file2,O", po::value<std::string>(&output_file2)->default_value(output_file2),
             "file for outputted match for reads of the second file to reads of the first one")
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
             "maximum distance value for trancated dist-graph construction")
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


        po::variables_map vm;
        store(po::command_line_parser(argc, argv).
              options(cmdline_options).run(), vm);
        notify(vm);

        if (config_file != "") {
            std::ifstream ifs(config_file.c_str());
            if (!ifs) {
                cout << "can not open config file: " << config_file << "\n";
                return 0;
            } else {
                store(parse_config_file(ifs, config_file_options), vm);
                // reparse cmd line again for update config defaults
                store(po::command_line_parser(argc, argv).
                      options(cmdline_options).run(), vm);
                notify(vm);
            }
        }

        if (vm.count("help-hidden")) {
            cout << cmdline_options << std::endl;
            return 0;
        }

        if (vm.count("help") || !vm.count("input-file1") || !vm.count("input-file2")) { // TODO Process required arguments by the proper way
            cout << visible << "\n";
            return 0;
        }

        if (vm.count("version")) {
            cout << "<Some cool name> version 0.1" << vm.count("version") << std::endl;
            return 0;
        }

        cout << "Input files are: "
            << vm["input-file1"].as<std::string>() << ", "<< vm["input-file2"].as<std::string>() << "\n";

        cout << "K = " << K << endl;
        cout << "tau = " << tau << endl;
    } catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
};

int main(int argc, char **argv) {
    cout << "Command line: " << join_cmd_line(argc, argv) << std::endl;


