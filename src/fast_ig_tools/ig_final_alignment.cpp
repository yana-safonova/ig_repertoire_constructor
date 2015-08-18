#define SEQAN_HAS_ZLIB 1 // TODO Set it using cmake
#define SEQAN_HAS_BZLIB 1


#include <fstream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <unordered_map>


#include <iostream>
using std::cout;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/format.hpp>
using bformat = boost::format;

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::SeqFileOut;
using seqan::CharString;

#include "ig_final_alignment.hpp"


int main(int argc, char **argv) {
  auto start_time = std::chrono::high_resolution_clock::now();

  int nthreads = 4;
  std::string input_file = "input.fa";
  std::string dense_subgraphs_file = "dense_subgraphs.txt";
  std::string output_file = "output.fa";
  std::string rcm_file = "repertoire.rcm";
  bool use_hamming_alignment = false;


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
      ("input-file,i", po::value<std::string>(&input_file),
       "name of the input file (FASTA|FASTQ)")
      ("dense-sgraphs,d", po::value<std::string>(&dense_subgraphs_file)->default_value(dense_subgraphs_file),
       "dense subgraphs file")
      ("output-file,o", po::value<std::string>(&output_file)->default_value(output_file),
       "output file for final repertoire")
      ("rcm-file,R", po::value<std::string>(&rcm_file)->default_value(rcm_file),
       "output RCM-file; empty string means no RCM-file producing")
      ("hamming,H",
       "use Hamming-based (position-wise) multiple alignment instead of seqan's one")
      ;

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
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

    if (vm.count("help") || !vm.count("input-file")) { // TODO Process required arguments by the proper way
      cout << visible << "\n";
      return 0;
    }

    if (vm.count("version")) {
      cout << "<Some cool name> version 0.1" << vm.count("version") << std::endl;
      return 0;
    }

    if (vm.count("hamming")) {
      use_hamming_alignment = true;
    }

    cout << "Input file is: "
      << input_file << "\n";
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  std::vector<Dna5String> input_reads;
  std::vector<CharString> input_ids; // Really, they are useless =)
  SeqFileIn seqFileIn_input(input_file.c_str());

  cout << "Reading data..." << std::endl;
  readRecords(input_ids, input_reads, seqFileIn_input);

  std::vector<size_t> component_indices;
  component_indices.reserve(input_reads.size());

  {
    std::ifstream dense_subgraphs_fh(dense_subgraphs_file.c_str());
    size_t index;
    while (dense_subgraphs_fh >> index) {
      component_indices.push_back(index);
    }
  }

  cout << bformat("Reads: %d; component_indices: %d") % input_reads.size() % component_indices.size() << std::endl;
  assert(component_indices.size() == input_reads.size());

  if (rcm_file != "") {
    cout << bformat("Saving RCM-file to %s...") % rcm_file << std::endl;
    std::ofstream rcm_fh(rcm_file);
    for (size_t i = 0; i < component_indices.size(); ++i) {
      rcm_fh << input_ids[i] << "\t" << component_indices[i] << "\n";
    }
  }

  size_t max_index = *std::max_element(component_indices.cbegin(), component_indices.cend());
  cout << bformat("The number of components: %d") % (max_index + 1) << std::endl;

  std::vector<std::vector<size_t>> component2id(max_index + 1);
  for (size_t i = 0; i < component_indices.size(); ++i) {
    component2id[component_indices[i]].push_back(i);
  }

  size_t max_component_size = 0;
  for (const auto &_ : component2id) {
    max_component_size = std::max(max_component_size, _.size());
  }
  cout << bformat("Maximum component size: %d") % max_component_size << std::endl;

  std::vector<Dna5String> consensuses(component2id.size());

  omp_set_num_threads(nthreads);
  cout << bformat("Consensus computation (using %d threads)...") % nthreads << std::endl;

  SEQAN_OMP_PRAGMA(parallel for)  // becomes: #pragma omp parallel for
  for (size_t comp = 0; comp < component2id.size(); ++comp) {
    if (component2id[comp].empty()) {
      continue;
    }

    if (use_hamming_alignment) {
      consensuses[comp] = consensus_hamming(input_reads, component2id[comp]);
    } else {
      consensuses[comp] = consensus(input_reads, component2id[comp]);
    }
  }

  cout << "Saving results" << std::endl;

  SeqFileOut seqFileOut_output(output_file.c_str());

  for (size_t comp = 0; comp < component2id.size(); ++comp) {
    if (component2id[comp].empty()) {
      continue;
    }

    size_t abundance = component2id[comp].size(); // TODO Use abundance from read ids

    bformat fmt("cluster___%d___size___%d");
    fmt % comp % abundance;
    std::string id = fmt.str();

    seqan::writeRecord(seqFileOut_output, id, consensuses[comp]);
  }

  cout << "Final repertoire has been saved to file " << output_file << std::endl;

  auto finish_time = std::chrono::high_resolution_clock::now();

  auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(finish_time - start_time).count();
  cout << bformat("Elapsed time: %0.3fs") % (double(elapsed_time) / 1000.) << std::endl;
  return 0;
}
