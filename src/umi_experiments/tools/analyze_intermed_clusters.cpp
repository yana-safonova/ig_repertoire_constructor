#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <segfault_handler.hpp>
#include <seqan/seq_io.h>
#include "utils.hpp"
#include "umi_utils.hpp"

namespace {
    struct Params {
        std::string reads_path;
        std::string cluster_centers_path;
        std::string rcm_path;
        std::string clusters_output_path;
        size_t cluster_size_threshold;
    };

    bool read_args(int argc, char **argv, Params& params) {
        namespace po = boost::program_options;
        po::options_description cmdl_options("Is this needed?");
        cmdl_options.add_options()
                ("help,h", "print help message")
                ("reads,r", po::value<std::string>(&params.reads_path)->required(), "input file with reads")
                ("clusters,c", po::value<std::string>(&params.cluster_centers_path)->required(), "file with reported clusters")
                ("rcm,m", po::value<std::string>(&params.rcm_path)->required(), "file with read-cluster map")
                ("cluster_reads,o", po::value<std::string>(&params.clusters_output_path)->required(), "directory to put large clusters to")
                ("size_threshold,s", po::value<size_t>(&params.cluster_size_threshold)->required(), "min cluster size to save")
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
        std::vector<seqan::CharString> input_ids;
        std::vector<seqan::Dna5String> input_reads;
        std::vector<seqan::CharString> cluster_ids;
        std::vector<seqan::Dna5String> cluster_reads;
        std::unordered_map<seqan::Dna5String, size_t> rcm;
    };

    Input read_everything(const Params& params) {
        Input input;

        INFO("Reading reads");
        seqan::SeqFileIn reads_file(params.reads_path.c_str());
        readRecords(input.input_ids, input.input_reads, reads_file);
        INFO(input.input_ids.size() << " records read");

        INFO("Reading clusters");
        seqan::SeqFileIn clusters_file(params.cluster_centers_path.c_str());
        readRecords(input.cluster_ids, input.cluster_reads, clusters_file);
        INFO(input.cluster_ids.size() << " clusters read");

        INFO("Reading rcm");
        std::ifstream rcm(params.rcm_path);
        while (!rcm.eof()) {
            std::string s;
            std::getline(rcm, s);
            if (s.empty()) continue;
            const std::vector<string>& tokens = split(s, '\t');
            VERIFY(tokens.size() == 2);
            const auto& id = seqan::CharString(tokens[0]);
            VERIFY(input.rcm.count(id) == 0);
            input.rcm[id] = std::stoull(tokens[1]);
        }

        return input;
    }
}

void save_clusters(const std::unordered_map<size_t, std::unordered_set<size_t>>& clusters, const std::vector<seqan::CharString>& cluster_ids,
                   const std::vector<seqan::CharString>& input_ids, const std::vector<seqan::Dna5String>& input_reads,
                   size_t threshold, const std::string& output_dir) {
    namespace fs = boost::filesystem;
    const fs::path output_dir_path(output_dir);
    if (fs::exists(output_dir_path)) {
        fs::remove_all(output_dir_path);
    }
    fs::create_directory(output_dir_path);

    for (const auto& cluster : clusters) {
        const auto& read_idxs = cluster.second;
        if (read_idxs.size() < threshold) continue;
        const std::string& cluster_file_name = seqan_string_to_string(cluster_ids[cluster.first]) + ".fasta";
        fs::path cluster_path(output_dir_path);
        cluster_path /= cluster_file_name;
        std::vector<seqan::CharString> read_ids;
        std::vector<seqan::Dna5String> reads;
        for (size_t idx : read_idxs) {
            read_ids.push_back(input_ids[idx]);
            reads.push_back(input_reads[idx]);
        }
        seqan::SeqFileOut cluster_file(cluster_path.string().c_str());
        writeRecords(cluster_file, read_ids, reads);
    }
}

int main(int argc, char **argv) {
    segfault_handler sh;
    create_console_logger();

    Params params;
    if (!read_args(argc, argv, params)) {
        return 0;
    }

    const auto& input = read_everything(params);

    std::unordered_map<seqan::CharString, size_t> read_id_to_idx;
    for (size_t i = 0; i < input.input_ids.size(); i ++) {
        read_id_to_idx[input.input_ids[i]] = i;
    }

    std::unordered_map<size_t, std::unordered_set<size_t>> clusters;
    for (const auto& entry : input.rcm) {
        size_t read_idx = read_id_to_idx[entry.first];
        clusters[entry.second].insert(read_idx);
    }

    INFO("Saving clusters of size at least " << params.cluster_size_threshold);
    save_clusters(clusters, input.cluster_ids, input.input_ids, input.input_reads, params.cluster_size_threshold, params.clusters_output_path);
    INFO("Done saving");
}