#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <segfault_handler.hpp>
#include <seqan/seq_io.h>
#include <clusterer.hpp>
#include "utils.hpp"
#include "umi_utils.hpp"

namespace {
    struct Params {
        std::string repertoire_path;
        std::string rcm_path;
        std::string reference_path;
        std::string reads_path;
    };

    bool read_args(int argc, char **argv, Params& params) {
        namespace po = boost::program_options;
        po::options_description cmdl_options("Is this needed?");
        cmdl_options.add_options()
                ("help,h", "print help message")
                ("repertoire,r", po::value<std::string>(&params.repertoire_path)->required(), "file with constructed repertoire")
                ("rcm,m", po::value<std::string>(&params.rcm_path)->required(), "file with read-cluster map for constructed repertoire")
                ("reference,f", po::value<std::string>(&params.reference_path)->required(), "file with reference repertoire")
                ("reads,s", po::value<std::string>(&params.reads_path)->required(), "file with original reads")
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
        std::vector<seqan::CharString> repertoire_ids;
        std::vector<seqan::Dna5String> repertoire_reads;
        std::unordered_map<size_t, std::vector<seqan::CharString>> rcm;
        std::vector<seqan::CharString> reference_ids;
        std::vector<seqan::Dna5String> reference;
        std::unordered_map<seqan::CharString, seqan::Dna5String> read_id_to_read;
    };

    Input read_everything(const Params& params) {
        Input input;

        {
            INFO("Reading constructed repertoire");
            seqan::SeqFileIn reads_file(params.repertoire_path.c_str());
            readRecords(input.repertoire_ids, input.repertoire_reads, reads_file);
            INFO(input.repertoire_ids.size() << " records read");
        }

        INFO("Reading rcm");
        std::ifstream rcm(params.rcm_path);
        while (!rcm.eof()) {
            std::string s;
            std::getline(rcm, s);
            if (s.empty()) continue;
            const std::vector<string>& tokens = split(s, '\t');
            VERIFY(tokens.size() == 2);
            const auto& id = seqan::CharString(tokens[0]);
            input.rcm[std::stoull(tokens[1])].push_back(id);
        }

        {
            INFO("Reading reference");
            seqan::SeqFileIn reference_file(params.reference_path.c_str());
            readRecords(input.reference_ids, input.reference, reference_file);
            INFO(input.reference.size() << " records read");
        }

        {
            INFO("Reading original reads");
            seqan::SeqFileIn original_reads_file(params.reads_path.c_str());
            std::vector<seqan::CharString> read_ids;
            std::vector<seqan::Dna5String> reads;
            readRecords(read_ids, reads, original_reads_file);
            for (size_t i = 0; i < reads.size(); i ++) {
                input.read_id_to_read[read_ids[i]] = reads[i];
            }
            INFO(reads.size() << " records read");
        }

        return input;
    }
}

//void save_clusters(const std::unordered_map<size_t, std::unordered_set<size_t>>& clusters, const std::vector<seqan::CharString>& cluster_ids,
//                   const std::vector<seqan::CharString>& input_ids, const std::vector<seqan::Dna5String>& input_reads,
//                   size_t threshold, const std::string& output_dir) {
//    namespace fs = boost::filesystem;
//    const fs::path output_dir_path(output_dir);
//    if (fs::exists(output_dir_path)) {
//        fs::remove_all(output_dir_path);
//    }
//    fs::create_directory(output_dir_path);
//
//    for (const auto& cluster : clusters) {
//        const auto& read_idxs = cluster.second;
//        if (read_idxs.size() < threshold) continue;
//        const std::string& cluster_file_name = seqan_string_to_string(cluster_ids[cluster.first]) + ".fasta";
//        fs::path cluster_path(output_dir_path);
//        cluster_path /= cluster_file_name;
//        std::vector<seqan::CharString> read_ids;
//        std::vector<seqan::Dna5String> reads;
//        for (size_t idx : read_idxs) {
//            read_ids.push_back(input_ids[idx]);
//            reads.push_back(input_reads[idx]);
//        }
//        seqan::SeqFileOut cluster_file(cluster_path.string().c_str());
//        writeRecords(cluster_file, read_ids, reads);
//    }
//}

void find_cluster(Input& input) {
    const auto dist0 = clusterer::ClusteringMode::bounded_edit_dist(0, 1);
    const auto dist = clusterer::ClusteringMode::bounded_edit_dist(2, 1, false);
    size_t bad_count = 0;
    for (size_t i = 0; i < input.repertoire_ids.size(); i ++) {
        bool found = false;
        for (size_t j = 0; j < input.reference.size(); j ++) {
            if (input.repertoire_reads[i] == input.reference[j]) {
                found = true;
                break;
            }
        }
        if (found) continue;
        for (size_t j = 0; j < input.reference.size(); j ++) {
            if (dist(input.repertoire_reads[i], input.reference[j]) == 1) {
                bad_count ++;

                auto id = seqan_string_to_string(input.repertoire_ids[i]);
                size_t size = stoull(id.substr(id.rfind('_') + 1));

                std::string cut_cluster = id.substr(id.find('_') + 3);
                const size_t cluster_id = stoull(cut_cluster.substr(0, cut_cluster.find('_')));
                INFO("Considering cluster with id " << cluster_id << " (" << input.repertoire_ids[i] << ")");
                const auto cluster_read_ids = input.rcm[cluster_id];

                if (cluster_read_ids.size() != size) {
                    FATAL_ERROR("Sizes diverge.");
                }
                if (cluster_read_ids.size() > 10) continue;

                for (const auto& read_id : cluster_read_ids) {
                    std::cout << ">" << read_id << "\n" << input.read_id_to_read[read_id] << "\n";
                }

                std::cout << ">" << input.repertoire_ids[i] << "_Consensus:\n" << input.repertoire_reads[i] << "\n<"
                          << input.reference_ids[j] << "_Reference:\n" << input.reference[j] << "\n\n\n";
                break;
//                return;
            }
        }
    }
    INFO("Total " << bad_count << " bad clusters found");
}

int main(int argc, char **argv) {
    segfault_handler sh;
    create_console_logger();

    Params params;
    if (!read_args(argc, argv, params)) {
        return 0;
    }

    auto input = read_everything(params);

    find_cluster(input);

}