#include <string>
#include "../umi_experiments/utils.hpp"
#include "../umi_experiments/umi_utils.hpp"
#include <segfault_handler.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <logger/log_writers.hpp>

namespace {

    struct Params {
        std::string input_file;
    };

    bool read_args(int argc, const char *const *argv, Params& params) {
        namespace po = boost::program_options;
        po::options_description cmdl_options;
        cmdl_options.add_options()
                ("help,h", "print help message")
                ("input,i", po::value<std::string>(&params.input_file)->required(), "input file with reads");
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(cmdl_options).run(), vm);
        if (vm.count("help") || argc == 1) {
            std::cout << cmdl_options << std::endl;
            return false;
        }
        po::notify(vm);
        return true;
    }

    void read_input(const Params& params, std::vector<seqan::CharString>& input_ids, std::vector<seqan::Dna5String>& input_reads) {
        INFO("Reading records from " << params.input_file);
        seqan::SeqFileIn reads_file(params.input_file.c_str());
        readRecords(input_ids, input_reads, reads_file);
        INFO(input_ids.size() << " records read");
    }
}

void report_umi_group_sizes(const std::unordered_map<Umi, std::vector<size_t>>& umi_to_reads) {
    std::map<size_t, size_t> size_to_cnt;
    for (const auto& entry : umi_to_reads) {
        size_to_cnt[entry.second.size()] ++;
    }
    for (const auto& entry : size_to_cnt) {
        TRACE("Size " << entry.first << " -> total " << entry.second);
    }
}

void report_umi_length_distribution(const std::unordered_map<Umi, std::vector<size_t>>& umi_to_reads) {
    std::map<size_t, size_t> umi_len_to_cnt;
    for (const auto& entry : umi_to_reads) {
        umi_len_to_cnt[length(entry.first.GetString())] ++;
    }
    for (const auto& entry : umi_len_to_cnt) {
        TRACE("Barcode length " << entry.first << " -> total " << entry.second);
    }
}

seqan::Dna5String compress_id(const seqan::Dna5String& id) {
    std::string id_str = seqan_string_to_string(id);
    std::string compressed = "";
    char prev = 0;
    for (char cur : id_str) {
        if (cur != prev) {
            compressed += cur;
        }
        prev = cur;
    }
    if (compressed.size() == 0) {
        ERROR("'" << id << "' -> '" << compressed << "'");
    }
    if (length(seqan::Dna5String(compressed)) == 0) {
        ERROR("'" << id << "' -> '" << compressed << "'");
    }
    return compressed;
}

int main(int argc, const char* const* argv) {
    segfault_handler sh;
    create_console_logger(logging::L_TRACE);

    Params params;
    if (!read_args(argc, argv, params)) {
        return 0;
    }

    std::vector<seqan::CharString> input_ids;
    std::vector<seqan::Dna5String> input_reads;
    read_input(params, input_ids, input_reads);

    std::vector<seqan::Dna5String> umis;
    extract_barcodes_from_read_ids(input_ids, umis);
    INFO("Total " << umis.size() << " barcodes extracted.")

    std::unordered_map<Umi, std::vector<size_t>> umi_to_reads;
    group_reads_by_umi(umis, umi_to_reads);

    report_umi_group_sizes(umi_to_reads);
    report_umi_length_distribution(umi_to_reads);

/*
    std::map<size_t, size_t> start_to_cnt;
    std::map<size_t, size_t> end_to_cnt;
    for (size_t i = 0; i < input_ids.size(); i ++) {
        const std::string s = seqan_string_to_string(input_reads[i]);
        start_to_cnt[s.find("TGAGCGGAAC")] ++;
        end_to_cnt[s.find("GGGAATTCTCACAGGAGACG")] ++;
    }
    for (const auto& entry : start_to_cnt) {
        INFO("Start pos " << entry.first << " -> total " << entry.second);
    }
    for (const auto& entry : end_to_cnt) {
        INFO("End pos " << entry.first << " -> total " << entry.second);
    }
*/

    std::vector<seqan::Dna5String> compressed_umis;
    compressed_umis.reserve(umis.size());
    for (const auto& umi : umis) {
        compressed_umis.push_back(compress_id(umi));
    }
    std::unordered_map<Umi, std::vector<size_t>> compressed_umi_to_reads;
    group_reads_by_umi(compressed_umis, compressed_umi_to_reads);

    report_umi_group_sizes(compressed_umi_to_reads);
    report_umi_length_distribution(compressed_umi_to_reads);

    return 0;
}
