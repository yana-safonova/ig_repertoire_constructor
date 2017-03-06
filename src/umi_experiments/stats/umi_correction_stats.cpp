#include "string"
#include "../../ig_tools/utils/string_tools.hpp"
#include <umi_utils.hpp>
#include <logger/log_writers.hpp>
#include <segfault_handler.hpp>
#include <boost/program_options.hpp>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <utils/io.hpp>

size_t get_simulated_read_id(const seqan::CharString& id) {
    std::string s = seqan_string_to_string(id);
    const std::string original_prefix = "original_";
    if (s.find(original_prefix) != std::string::npos) {
        s = s.substr(original_prefix.length());
    }
    return std::stoull(s.substr(0, s.find('_')));
}

namespace {
    void create_console_logger() {
        using namespace logging;
        logger* lg = create_logger("");
        lg->add_writer(std::make_shared<console_writer>());
        attach_logger(lg);
    }

    struct Params {
        std::string amplified_reads_path;
        std::string amplification_rcm_path;
        std::string ideal_rcm_path;
        std::string igrec_rep_rcm_path;
    };

    bool read_args(int argc, const char* const* argv, Params& params) {
        namespace po = boost::program_options;
        po::options_description cmdl_options;
        cmdl_options.add_options()
                ("reads,r", po::value<std::string>(&params.amplified_reads_path)->required(), "path to a file with amplified reads")
                ("amplify-rcm,a", po::value<std::string>(&params.amplification_rcm_path)->required(), "path to an amplified read to original reads map")
                ("ideal-rcm,i", po::value<std::string>(&params.ideal_rcm_path)->required(), "path to an amplified read to original equal read cluster map")
                ("igrec-rcm,c", po::value<std::string>(&params.igrec_rep_rcm_path)->required(), "path to an amplified read to igrec constructed repertoire sequence map")
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
        Input(const std::vector<seqan::CharString>& amplified_read_ids_,
              const std::unordered_map<seqan::CharString, size_t>& amplified_to_original_,
              const std::unordered_map<seqan::CharString, size_t>& amplified_to_compressed_,
              const std::unordered_map<seqan::CharString, size_t>& igrec_rcm_)
                : amplified_read_ids(amplified_read_ids_),
                  amplified_to_original(amplified_to_original_),
                  amplified_to_compressed(amplified_to_compressed_),
                  igrec_rcm(igrec_rcm_) {}

        std::vector<seqan::CharString> amplified_read_ids;
        std::unordered_map<seqan::CharString, size_t> amplified_to_original;
        std::unordered_map<seqan::CharString, size_t> amplified_to_compressed;
        std::unordered_map<seqan::CharString, size_t> igrec_rcm;
    };

    Input read_input(const Params& params) {
        std::vector<seqan::CharString> amplified_read_ids_shuffled;
        {
            std::vector<seqan::Dna5String> amplified_reads;
            INFO("Reading records from " << params.amplified_reads_path);
            seqan::SeqFileIn reads_file(params.amplified_reads_path.c_str());
            readRecords(amplified_read_ids_shuffled, amplified_reads, reads_file);
//            amplified_read_ids_str = std::vector<std::string>(amplified_read_ids_shuffled.size());
//            std::transform(amplified_read_ids_shuffled.begin(), amplified_read_ids_shuffled.end(), amplified_read_ids_str.begin(), [](seqan::CharString s) -> std::string { return seqan_string_to_string(s); } );
            INFO(amplified_read_ids_shuffled.size() << " records read");
        }
        std::vector<seqan::CharString> amplified_read_ids(amplified_read_ids_shuffled.size());
        for (size_t i = 0; i < amplified_read_ids_shuffled.size(); i ++) {
            size_t id = get_simulated_read_id(amplified_read_ids_shuffled[i]);
            amplified_read_ids[id] = amplified_read_ids_shuffled[i];
        }
        const auto& amplified_to_original = read_rcm_file(params.amplification_rcm_path);
        const auto& amplified_to_compressed = read_rcm_file(params.ideal_rcm_path);
        const auto& igrec_rcm = read_rcm_file(params.igrec_rep_rcm_path);
        
        return Input(amplified_read_ids, amplified_to_original, amplified_to_compressed, igrec_rcm);
    }
}

int main(int argc, const char* const* argv) {
    segfault_handler sh;
    perf_counter pc;
    create_console_logger();

    Params params;
    if (!read_args(argc, argv, params)) return 0;
    Input input = read_input(params);

    std::vector<seqan::Dna5String> amplified_umis;
    extract_barcodes_from_read_ids(input.amplified_read_ids, amplified_umis);
    std::vector<size_t> reads_with_corrupted_barcode;
    for (size_t i = 0; i < input.amplified_read_ids.size(); i ++) {
        if (amplified_umis[i] != amplified_umis[input.amplified_to_original[input.amplified_read_ids[i]]]) {
            reads_with_corrupted_barcode.push_back(i);
        }
    }
    INFO(reads_with_corrupted_barcode.size() << " reads with error(s) in barcode generated");

    std::unordered_map<size_t, std::vector<seqan::CharString>> ideal_clusters;
    for (const auto& entry : input.amplified_to_compressed) {
        ideal_clusters[entry.second].push_back(entry.first);
    }
    INFO("Total ideal clusters: " << ideal_clusters.size());
    std::unordered_map<size_t, size_t> ideal_to_igrc_cluster;
    for (const auto& entry : ideal_clusters) {
        std::unordered_map<size_t, size_t> igrc_cluster_to_size;
        for (const auto& id : entry.second) {
            igrc_cluster_to_size[input.igrec_rcm[id]] ++;
        }
        size_t max = 0;
        for (const auto& igrc_cluster_and_size : igrc_cluster_to_size) {
            if (igrc_cluster_and_size.second > max) {
                max = igrc_cluster_and_size.second;
                ideal_to_igrc_cluster[entry.first] = igrc_cluster_and_size.first;
            }
        }
    }

    size_t corrected = 0;
    for (size_t read_idx : reads_with_corrupted_barcode) {
        const seqan::CharString read = input.amplified_read_ids[read_idx];
        if (input.igrec_rcm[read] == ideal_to_igrc_cluster[input.amplified_to_compressed[read]]) {
            corrected ++;
        }
    }
    INFO(corrected << " of them have barcode corrected.");
}
