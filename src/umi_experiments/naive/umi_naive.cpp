#include <utils.hpp>
#include <boost/program_options.hpp>
#include <umi_utils.hpp>
#include <clusterer.hpp>
#include <segfault_handler.hpp>

namespace {
    struct Params {
        std::string reads_path;
        std::string umi_uncompressed_path;
        std::string umi_compressed_path;
        std::string umi_graph_path;
        std::string output_dir;
        bool save_clusters;
    };

    bool read_args(int argc, const char* const* argv, Params& params) {
        namespace po = boost::program_options;
        po::options_description cmdl_options("Is this needed?");
        cmdl_options.add_options()
                ("reads,r", po::value<std::string>(&params.reads_path)->required(), "input file with reads")
                ("output,o", po::value<std::string>(&params.output_dir)->default_value(""), "output directory path")
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
        Input(std::vector<seqan::CharString>& input_ids_, std::vector<seqan::Dna5String>& input_reads_)
                : input_ids(input_ids_), input_reads(input_reads_)
        {}

        std::vector<seqan::CharString> input_ids;
        std::vector<seqan::Dna5String> input_reads;
    };

    Input read_everything(const Params& params) {
        std::vector<seqan::CharString> input_ids;
        std::vector<seqan::Dna5String> input_reads;
        INFO("Reading records from " << params.reads_path);
        seqan::SeqFileIn reads_file(params.reads_path.c_str());
        readRecords(input_ids, input_reads, reads_file);
        INFO(input_ids.size() << " records read");

        return Input(input_ids, input_reads);
    }
}

std::vector<string> extract_umis(const std::vector<seqan::CharString> read_ids) {
    std::vector<seqan::Dna5String> umis_seqan;
    std::vector<seqan::DnaQString> umi_quals;
    extract_barcodes_from_read_ids(read_ids, umis_seqan, umi_quals);
    std::vector<string> umis(read_ids.size());
    std::transform(umis_seqan.cbegin(), umis_seqan.cend(), umis.begin(), seqan_string_to_string<seqan::Dna5String>);
    return umis;
}

seqan::Dna5String calculate_consensus(const std::vector<seqan::Dna5String>& reads, const std::vector<size_t>& idx_list) {
    VERIFY_MSG(!idx_list.empty(), "Calculating consensus by empty index list");
    seqan::Dna5String consensus(reads[idx_list[0]]);
    for (size_t idx : idx_list) {
        if (length(reads[idx]) < length(consensus)) {
            consensus = seqan::Dna5String(reads[idx]);
        }
    }

    auto GetIndexNucl = [](const seqan::Dna5 &nucl) {
        return static_cast<size_t>(seqan::ordValue(nucl));
    };

    for (size_t pos = 0; pos < length(consensus); pos ++) {
        std::vector<size_t> nt_counts(4);
        for (size_t idx : idx_list) {
            nt_counts[GetIndexNucl(reads[idx][pos])] ++;
        }
        for (size_t i = 0; i < 4; i ++) {
            if (nt_counts[i] > nt_counts[GetIndexNucl(consensus[pos])]) {
                consensus[pos] = i;
            }
        }
    }
    return consensus;
}

int main(int argc, const char* const* argv) {
    segfault_handler sh;
    create_console_logger();

    Params params;
    if (!read_args(argc, argv, params)) {
        return 0;
    }

    const auto& input = read_everything(params);
    const std::vector<string>& umis = extract_umis(input.input_ids);
    
    std::unordered_map<std::string, std::vector<size_t>> umi_to_idx_list;
    for (size_t i = 0; i < input.input_ids.size(); i ++) {
        auto& idx_list = umi_to_idx_list[umis[i]];
        idx_list.push_back(i);
    }

    std::unordered_map<std::string, seqan::Dna5String> umi_to_consensus;
    size_t min_consensus_length = std::numeric_limits<size_t>::max();
    size_t max_consensus_length = 0;
    for (auto& entry : umi_to_idx_list) {
        const auto& idx_list = entry.second;
        const auto& consensus = calculate_consensus(input.input_reads, idx_list);
        std::vector<size_t> filtered_idx_list;
        const auto& dist = clusterer::ClusteringMode::bounded_hamming_dist(10);
        std::remove_copy_if(idx_list.begin(), idx_list.end(), back_inserter(filtered_idx_list), [&consensus, &dist](seqan::Dna5String read) {
            return dist(consensus, read) > 10;
        });
        entry.second = filtered_idx_list;
        const auto& precise_consensus = calculate_consensus(input.input_reads, filtered_idx_list);
        umi_to_consensus[entry.first] = precise_consensus;
        min_consensus_length = min(min_consensus_length, length(precise_consensus));
        max_consensus_length = max(max_consensus_length, length(precise_consensus));
    }

    std::unordered_map<std::string, std::vector<std::string>> consensus_prefix_to_umi_list;
    for (const auto& entry : umi_to_consensus) {
        const std::string& prefix = seqan_string_to_string(entry.second).substr(0, min_consensus_length);
        consensus_prefix_to_umi_list[prefix].push_back(entry.first);
    }

    std::vector<size_t> cluster_id(input.input_reads.size(), std::numeric_limits<size_t>::max());
    std::vector<seqan::Dna5String> repertoire;
    std::vector<seqan::CharString> repertoire_ids;
    for (const auto& entry : consensus_prefix_to_umi_list) {
        std::vector<size_t> idx_list;
        for (const auto& umi : entry.second) {
            const vector<size_t>& list = umi_to_idx_list[umi];
            idx_list.insert(idx_list.end(), list.begin(), list.end());
        }
        const seqan::Dna5String& consensus = calculate_consensus(input.input_reads, idx_list);
        size_t sequence_id = repertoire.size();
        for (size_t idx : idx_list) {
            cluster_id[idx] = sequence_id;
        }
        repertoire.push_back(consensus);
        repertoire_ids.push_back(seqan::CharString("cluster___" + to_string(sequence_id) + "___size___" + to_string(idx_list.size())));
    }
    write_seqan_records(boost::filesystem::path(params.output_dir).append("naive_repertoire.fa"), repertoire_ids, repertoire);
    std::ofstream rcm_file(boost::filesystem::path(params.output_dir).append("naive_repertoire.rcm").string());
    for (size_t idx = 0; idx < input.input_reads.size(); idx ++) {
        rcm_file << input.input_ids[idx] << "\t" << (cluster_id[idx] < repertoire.size() ? to_string(cluster_id[idx]) : "") << std::endl;
    }

    return 0;
}
