#include <chrono>

#include <future>
using std::async;

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

#include <openmp_wrapper.h>

#include "ig_matcher.hpp"
#include "banded_half_smith_waterman.hpp"
#include "fast_ig_tools.hpp"


class BestScoreIndices {
    public:
        std::vector<size_t> indices;
        static const int INF = 1005000;
        int best_score;

        BestScoreIndices() : best_score{-INF} {
            omp_init_lock(&lock);
        }

        ~BestScoreIndices() {
            omp_destroy_lock(&lock);
        }

        void update(int score, size_t index) {
            omp_set_lock(&lock);
            update_nolock(score, index);
            omp_unset_lock(&lock);
        }

        void update_nolock(int score, size_t index) {
            if (score < best_score) {
                // Do nothing
            } else if (score == best_score) {
                indices.push_back(index);
            } else {
                best_score = score;
                indices = std::vector<size_t>( { index } );
            }
        }

    private:
        omp_lock_t lock;
};


// Format: in each line <best_score> id1, id2, ...
void write_best_hits(const std::vector<BestScoreIndices> &g, const std::string &filename) {
    std::ofstream out(filename);

    for (const auto &_ : g) {
        if (_.indices.empty()) {
            out << "\n";
            continue;
        }

        out << _.best_score;
        for (size_t i : _.indices) {
            out << " " << i;
        }
        out << "\n";
    }
}


template<typename T, typename Tf>
void bestScorePairing(const std::vector<T> &input_reads1,
                      const std::vector<T> &input_reads2,
                      const KmerIndex &kmer2reads1,
                      const KmerIndex &kmer2reads2,
                      const Tf &score_fun,
                      int tau,
                      int K,
                      Strategy strategy,
                      const std::string &fname1,
                      const std::string &fname2) {
    std::vector<BestScoreIndices> g1(input_reads1.size());
    std::vector<BestScoreIndices> g2(input_reads2.size());

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 8))
        for (size_t j = 0; j < input_reads1.size(); ++j) {
            auto cand = find_candidates(input_reads1[j], kmer2reads2, input_reads2.size(), tau, K, strategy);

            for (size_t i : cand) {
                int score = score_fun(input_reads1[j], input_reads2[i]);
                g1[j].update_nolock(score, i);
                g2[i].update(score, j); // With lock !!!
            }
        }

    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic, 8))
        for (size_t j = 0; j < input_reads2.size(); ++j) {
            auto cand = find_candidates(input_reads2[j], kmer2reads1, input_reads1.size(), tau, K, strategy);

            for (size_t i : cand) {
                int score = score_fun(input_reads2[j], input_reads1[i]);
                g2[j].update_nolock(score, i);
                g1[i].update(score, j); // With lock !!!
            }
        }

    SEQAN_OMP_PRAGMA(parallel for schedule(guided, 8))
        for (size_t i = 0; i < g1.size(); ++i) {
            remove_duplicates(g1[i].indices);
        }

    SEQAN_OMP_PRAGMA(parallel for schedule(guided, 8))
        for (size_t i = 0; i < g2.size(); ++i) {
            remove_duplicates(g2[i].indices);
        }

    write_best_hits(g1, fname1);
    write_best_hits(g2, fname2);
}


int main(int argc, char **argv) {
    auto start_time = std::chrono::high_resolution_clock::now();

    cout << "Command line: " << join_cmd_line(argc, argv) << std::endl;

    int K = 36; // anchor length
    int tau = 3;
    int nthreads = 4;
    std::string input_file1 = "input1.fa";
    std::string input_file2 = "input2.fa";
    std::string output_file1 = "output1.match";
    std::string output_file2 = "output2.match";
    int strategy_int = Strategy::TRIPLE;

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

    SeqFileIn seqFileIn_input1(input_file1.c_str());
    SeqFileIn seqFileIn_input2(input_file2.c_str());
    std::vector<CharString> input_ids1, input_ids2; // Really, they are useless =)
    std::vector<Dna5String> input_reads1, input_reads2;

    cout << "Reading data..." << std::endl;
    readRecords(input_ids1, input_reads1, seqFileIn_input1);
    readRecords(input_ids2, input_reads2, seqFileIn_input2);
    cout << bformat("Reads: %d\n") % length(input_reads1);
    cout << bformat("Reads: %d\n") % length(input_reads2);

    cout << "Reads' length checking..." << std::endl;
    size_t required_read_length = (strategy_int != 0) ? (K * (tau + strategy_int)) : 0;

    size_t discarded_reads1 = 0;
    for (const auto &read : input_reads1) {
        discarded_reads1 += length(read) < required_read_length;
    }
    size_t discarded_reads2 = 0;
    for (const auto &read : input_reads2) {
        discarded_reads2 += length(read) < required_read_length;
    }

    if (discarded_reads1 || discarded_reads2) {
        cout << bformat("Discarded reads %d/%d") % discarded_reads1 % discarded_reads2 << std::endl;
    }

    cout << "K-mer index construction (||)..." << std::endl;
    // auto kmer2reads1 = kmerIndexConstruction(input_reads1, K);
    // auto kmer2reads2 = kmerIndexConstruction(input_reads2, K);
    auto kmer2reads1_fut = async(std::launch::async, [&](){ return kmerIndexConstruction(input_reads1, K); });
    auto kmer2reads2_fut = async(std::launch::async, [&](){ return kmerIndexConstruction(input_reads2, K); });

    auto kmer2reads1 = kmer2reads1_fut.get();
    auto kmer2reads2 = kmer2reads2_fut.get();

    omp_set_num_threads(nthreads);
    cout << bformat("Matching (using %d threads)...") % nthreads << std::endl;
    Strategy strategy = Strategy(strategy_int);
    cout << toCString(strategy) << std::endl;

    /*
    auto score_fun = [tau](const Dna5String &s1, const Dna5String &s2) -> int {
        auto lizard_tail = [](int l) -> int { return ((l > 0) ? -1 : 0); };
        return half_sw_banded(s1, s2, 1, -1, -1, lizard_tail, tau);
    };
    */

    auto score_fun = [tau](const Dna5String &s1, const Dna5String &s2) -> int {
        auto lizard_tail = [](int) -> int { return 0; };
        return half_sw_banded(s1, s2, 0, -1, -1, lizard_tail, tau);
    };

    bestScorePairing(input_reads1, input_reads2,
                     kmer2reads1, kmer2reads2,
                     score_fun,
                     tau, K,
                     strategy,
                     output_file1, output_file2);

    cout << "Graph has been saved to files " << output_file1 << ", " << output_file2 << std::endl;

    auto finish_time = std::chrono::high_resolution_clock::now();

    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(finish_time - start_time).count();
    cout << bformat("Elapsed time: %0.3fs") % (double(elapsed_time) / 1000.) << std::endl;
    return 0;
}

// vim: ts=4:sw=4
