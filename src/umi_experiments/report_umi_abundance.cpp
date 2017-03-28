#include <string>
#include <unordered_map>
#include <unordered_set>

#include <segfault_handler.hpp>
#include <logger/log_writers.hpp>
#include <boost/program_options.hpp>
#include <seqan/file.h>
#include <utils/io.hpp>
#include "../ig_tools/utils/string_tools.hpp"
#include "umi_utils.hpp"

namespace {
    void create_console_logger() {
        using namespace logging;
        logger* lg = create_logger("");
        lg->add_writer(std::make_shared<console_writer>());
        attach_logger(lg);
    }

    struct Params {
        std::string final_rcm_path;
        std::string intermediate_rcm_path;
        std::string repertoire_path;
        std::string output_path;
    };

    bool read_args(int argc, const char* const* argv, Params& params) {
        namespace po = boost::program_options;
        po::options_description cmdl_options;
        cmdl_options.add_options()
                ("help,h", "print help message")
                ("final-rcm,f", po::value<std::string>(&params.final_rcm_path)->required(), "path to a file with final read-cluster map")
                ("inter-rcm,i", po::value<std::string>(&params.intermediate_rcm_path)->required(), "path to a file with intermediate read-cluster map")
                ("repertoire,r", po::value<std::string>(&params.repertoire_path)->required(), "path to a file with final repertoire")
                ("output,o", po::value<std::string>(&params.output_path)->required(), "path to an output file")
                ;
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(cmdl_options).run(), vm);
        if (vm.count("help") || argc == 1) {
            std::cout << cmdl_options << std::endl;
            return false;
        }
        po::notify(vm);
        return true;
    }

    struct Input {
        std::unordered_map<seqan::CharString, size_t> final_rcm;
        std::unordered_map<seqan::CharString, size_t> intermediate_rcm;
        std::vector<seqan::CharString> read_ids;
        std::vector<seqan::Dna5String> read_seqs;

        Input(std::unordered_map<seqan::CharString, size_t>& final_rcm_, std::unordered_map<seqan::CharString, size_t>& intermediate_rcm_,
              std::vector<seqan::CharString> read_ids_, std::vector<seqan::Dna5String> read_seqs_) :
                final_rcm(final_rcm_), intermediate_rcm(intermediate_rcm_), read_ids(read_ids_), read_seqs(read_seqs_) {}
    };

    Input read_everything(const Params& params) {
        auto final_rcm = read_rcm_file(params.final_rcm_path);
        auto intermediate_rcm = read_rcm_file(params.intermediate_rcm_path);
        std::vector<seqan::CharString> read_ids;
        std::vector<seqan::Dna5String> read_seqs;
        read_seqan_records(params.repertoire_path, read_ids, read_seqs);
        return Input(final_rcm, intermediate_rcm, read_ids, read_seqs);
    }

    std::unordered_map<size_t, size_t> get_umi_abundances(std::unordered_map<seqan::CharString, size_t>& intermediate_rcm,
                                                          std::unordered_map<seqan::CharString, size_t>& final_rcm) {
        VERIFY(intermediate_rcm.size() == final_rcm.size());
        std::unordered_map<size_t, std::unordered_set<size_t>> final_cluster_to_umi_clusters;
        std::unordered_map<size_t, size_t> umi_cluster_to_final_cluster;
        for (auto& entry : final_rcm) {
            const auto read_id = entry.first;
            const auto final_cluster = entry.second;
            VERIFY(intermediate_rcm.count(read_id));
            const auto intermed_cluster = intermediate_rcm[read_id];
            if (umi_cluster_to_final_cluster.count(intermed_cluster)) {
                VERIFY(umi_cluster_to_final_cluster[intermed_cluster] == final_cluster);
            }
            umi_cluster_to_final_cluster[intermed_cluster] = final_cluster;
            final_cluster_to_umi_clusters[final_cluster].insert(intermed_cluster);
        }
        std::unordered_map<size_t, size_t> final_cluster_to_size;
        for (auto& entry : final_cluster_to_umi_clusters) {
            final_cluster_to_size[entry.first] = entry.second.size();
        }
        return final_cluster_to_size;
    }
}

int main(int argc, const char* const* argv) {
    segfault_handler sh;
    perf_counter pc;
    create_console_logger();

    Params params;
    if (!read_args(argc, argv, params)) {
        return 0;
    }
    auto input = read_everything(params);

    INFO("Updating read ids");
    std::unordered_map<size_t, size_t> final_cluster_umi_abundances = get_umi_abundances(input.intermediate_rcm, input.final_rcm);
    std::vector<seqan::CharString> new_read_ids(input.read_ids.size());
    std::vector<size_t> cluster_umi_abundances(new_read_ids.size());
    for (size_t i = 0; i < input.read_ids.size(); i ++) {
        const std::string id = seqan_string_to_string(input.read_ids[i]);
        std::vector<string> parts = split(id, "___");
        size_t cluster = std::stoull(parts[1]);
//        size_t size = std::stoull(parts[3]);
//        new_read_ids[i] = seqan::CharString("cluster_" + std::to_string(cluster) + "|UMIs_" + std::to_string(final_cluster_umi_abundances[cluster]) + "|reads_" + std::to_string(size));
        size_t umi_abundances = final_cluster_umi_abundances[cluster];
        new_read_ids[i] = seqan::CharString("cluster___" + std::to_string(cluster) + "___size___" + std::to_string(umi_abundances));
        cluster_umi_abundances[i] = umi_abundances;
    }

    std::vector<size_t> sort_permutation(input.read_ids.size());
    std::iota(sort_permutation.begin(), sort_permutation.end(), 0);
    std::sort(sort_permutation.begin(), sort_permutation.end(), [&cluster_umi_abundances] (size_t left, size_t right) { return cluster_umi_abundances[left] > cluster_umi_abundances[right]; });

    std::vector<seqan::CharString> sorted_read_ids(input.read_ids.size());
    std::vector<seqan::Dna5String> sorted_read_seqs(sorted_read_ids.size());
    for (size_t i = 0; i < sorted_read_ids.size(); i ++) {
        sorted_read_ids[i] = new_read_ids[sort_permutation[i]];
        sorted_read_seqs[i] = input.read_seqs[sort_permutation[i]];
    }

//    write_seqan_records(params.output_path, new_read_ids, input.read_seqs);
    write_seqan_records(params.output_path, sorted_read_ids, sorted_read_seqs);
}
