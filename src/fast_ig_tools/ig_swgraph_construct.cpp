#define SEQAN_HAS_ZLIB 1 // TODO Set it using cmake
#define SEQAN_HAS_BZLIB 1

#include <chrono>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/format.hpp>
using bformat = boost::format;

#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#include <seqan/seq_io.h>
using seqan::Dna5String;
using seqan::SeqFileIn;
using seqan::CharString;

#include "ig_matcher.hpp"
#include "banded_half_smith_waterman.hpp"


template<typename T, typename Tf>
Graph tauDistGraph(const std::vector<T> &input_reads,
                   const KmerIndex &kmer2reads,
                   const Tf &dist_fun,
                   int tau,
                   int K,
                   Strategy strategy) {
  Graph g(input_reads.size());

  SEQAN_OMP_PRAGMA(parallel for)  // becomes: #pragma omp parallel for
  for (size_t j = 0; j < input_reads.size(); ++j) {
    auto cand = find_candidates(input_reads[j], kmer2reads, input_reads.size(), tau, K, strategy);

    for (size_t i : cand) {
      int dist = dist_fun(input_reads[j], input_reads[i]);
      if (dist <= tau) {
        g[j].push_back( { i, dist } );
      }
    }
  }

  // Undirecting
  auto gg = g;
  for (size_t i = 0; i < gg.size(); ++i) {
    for (const auto &_ : gg[i]) {
      g[_.first].push_back( { i, _.second } );
    }
  }
  gg.clear(); // Free memory

  SEQAN_OMP_PRAGMA(parallel for)  // becomes: #pragma omp parallel for
  for (size_t j = 0; j < g.size(); ++j) {
    remove_duplicates(g[j]);
  }

  return g;
}


int main(int argc, char **argv) {
  auto start_time = std::chrono::high_resolution_clock::now();

  int K = 36; // anchor length
  int tau = 3;
  int nthreads = 4;
  std::string input_file = "cropped.fa";
  std::string output_filename = "output.graph";
  int strategy_int = Strategy::TRIPLE;
  int max_indels = 0;

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
       "name of an input file (FASTA|FASTQ)")
      ("output-file,o", po::value<std::string>(&output_filename)->default_value(output_filename),
       "file for outputted truncated dist-graph in METIS format")
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
      ("max-indels", po::value<int>(&max_indels)->default_value(max_indels),
       "maximum number of indels in Levenshtein distance")
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

    po::positional_options_description p;
    p.add("input-file", -1);
    // p.add("output-file", -1);

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

    cout << "K = " << K << endl;
    cout << "tau = " << tau << endl;
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }

  SeqFileIn seqFileIn_input(input_file.c_str());
  std::vector<CharString> input_ids;
  std::vector<Dna5String> input_reads;

  cout << "Reading data..." << std::endl;
  readRecords(input_ids, input_reads, seqFileIn_input);
  cout << bformat("Reads: %d\n") % length(input_reads);

  cout << "Reads' length checking..." << std::endl;
  size_t required_read_length = (strategy_int != 0) ? (K * (tau + strategy_int)) : 0;

  size_t discarded_reads = 0;
  for (const auto &read : input_reads) {
    discarded_reads += length(read) < required_read_length;
  }

  if (discarded_reads) {
    cout << bformat("Discarded reads %d") % discarded_reads << std::endl;
  }

  cout << "K-mer index construction..." << std::endl;
  auto kmer2reads = kmerIndexConstruction(input_reads, K);

  omp_set_num_threads(nthreads);
  cout << bformat("Truncated dist-graph construction (using %d threads)") % nthreads << std::endl;
  cout << "Candidates graph construcion..." << std::endl;

  Strategy strategy = Strategy(strategy_int);
  cout << toCString(strategy) << std::endl;

  auto dist_fun = [max_indels](const Dna5String &s1, const Dna5String &s2) -> int {
    auto lizard_tail = [](int l) -> int { return 0*l; };
    return -half_sw_banded(s1, s2, 0, -1, -1, lizard_tail, max_indels);
  };

  auto dist_graph = tauDistGraph(input_reads,
                                 kmer2reads,
                                 dist_fun,
                                 tau, K,
                                 strategy);

  // Output
  cout << "Saving graph..." << std::endl;
  write_metis_graph(dist_graph, output_filename);
  cout << "Graph has been saved to file " << output_filename << std::endl;

  auto finish_time = std::chrono::high_resolution_clock::now();

  auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(finish_time - start_time).count();
  cout << bformat("Elapsed time: %0.3fs") % (double(elapsed_time) / 1000.) << std::endl;
  return 0;
}
