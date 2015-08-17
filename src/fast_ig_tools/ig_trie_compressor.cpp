#define SEQAN_HAS_ZLIB 1 // TODO Set it using cmake
#define SEQAN_HAS_BZLIB 1


#include <vector>
#include <string>
#include <cassert>
#include <algorithm>
#include <chrono>

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::SeqFileOut;
using seqan::CharString;
using seqan::length;


#include <boost/format.hpp>
using bformat = boost::format;

#include <iostream>
using std::cout;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "ig_trie_compressor.hpp"


int main(int argc, char **argv) {
  auto start_time = std::chrono::high_resolution_clock::now();

  std::string input_file = "input.fa";
  std::string output_file = "output.fa";
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
      ("input-file,i", po::value<std::string>(&input_file),
       "name of the input file (FASTA|FASTQ)")
      ("output-file,o", po::value<std::string>(&output_file)->default_value(output_file),
       "name of the output file (FASTA|FASTQ)")
      ;

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
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

    po::positional_options_description p;
    p.add("input-file", -1);

    po::variables_map vm;
    store(po::command_line_parser(argc, argv).
          options(cmdline_options).positional(p).run(), vm);
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
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);
      }
    }

    if (vm.count("help-hidden")) {
      cout << cmdline_options << std::endl;
      return 0;
    }

    if (vm.count("help") || !vm.count("input-file")) { // TODO Process required arguments by the proper way
      cout << visible << "\n";
      return 0;
    }

    if (vm.count("version")) {
      cout << "<Some cool name> version 0.1" << vm.count("version") << std::endl;
      return 0;
    }

    if (vm.count("input-file")) {
      cout << "Input file is: "
        << vm["input-file"].as<std::string>() << "\n";
    }

    if (vm.count("output-file")) {
      cout << "Output file is: "
        << vm["output-file"].as<std::string>() << "\n";
    }
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  SeqFileIn seqFileIn_input(input_file.c_str());
  SeqFileOut seqFileOut_output(output_file.c_str());
  std::vector<CharString> input_ids;
  std::vector<Dna5String> input_reads;

  cout << "Reading data..." << std::endl;
  readRecords(input_ids, input_reads, seqFileIn_input);
  cout << bformat("Reads: %d\n") % length(input_reads);

  Trie<5> trie;

  cout << "Construction trie..." << std::endl;
  for (size_t i = 0; i < length(input_reads); ++i) {
    trie.add(input_reads[i], i);
  }

  cout << "Unique prefixes collecting..." << std::endl;
  auto result = trie.checkout(length(input_reads));

  cout << "Saving output..." << std::endl;
  size_t count = 0;
  for (const auto &_ : result) {
    size_t index = _.first;
    size_t abundance = _.second;

    count += abundance;

    std::string id = seqan::toCString(input_ids[index]);
    id += "_abundance:" + std::to_string(abundance);

    seqan::writeRecord(seqFileOut_output, id, input_reads[index]);
  }

  assert(count == length(input_reads));

  cout << bformat("Output reads: %d\n") % result.size();

  auto finish_time = std::chrono::high_resolution_clock::now();

  auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(finish_time - start_time).count();
  cout << bformat("Elapsed time: %0.3fs") % (double(elapsed_time) / 1000.) << std::endl;
  return 0;
}
