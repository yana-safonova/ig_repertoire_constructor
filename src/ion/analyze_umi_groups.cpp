#include <string>
#include "../umi_experiments/utils.hpp"
#include "../umi_experiments/umi_utils.hpp"
#include "ion_utils.hpp"
#include <segfault_handler.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <logger/log_writers.hpp>

namespace {

    struct Params {
        std::string input_file;
        std::string output_file;
    };

    bool read_args(int argc, const char *const *argv, Params& params) {
        namespace po = boost::program_options;
        po::options_description cmdl_options;
        cmdl_options.add_options()
                ("help,h", "print help message")
                ("input,i", po::value<std::string>(&params.input_file)->required(), "input file with reads")
                ("compressed,c", po::value<std::string>(&params.output_file)->required(), "output file with compressed barcodes");
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(cmdl_options).run(), vm);
        if (vm.count("help") || argc == 1) {
            std::cout << cmdl_options << std::endl;
            return false;
        }
        po::notify(vm);
        return true;
    }

    void read_input(const std::string& input_file,
                    std::vector<seqan::CharString>& ids,
                    std::vector<seqan::Dna5String>& reads) {
        INFO("Reading records from " << input_file);
        seqan::SeqFileIn reads_file(input_file.c_str());
        readRecords(ids, reads, reads_file);
        INFO(ids.size() << " records read");
    }

    void write_output(const std::string output_path,
                      const std::vector<seqan::CharString>& ids,
                      const std::vector<seqan::Dna5String>& reads) {
        seqan::SeqFileOut output_file(output_path.c_str());
        writeRecords(output_file, ids, reads);
        INFO(ids.size() << " sequences written to " << output_path << ".");
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
    read_input(params.input_file, input_ids, input_reads);

    std::vector<seqan::Dna5String> umis;
    extract_barcodes_from_read_ids(input_ids, umis);
    INFO("Total " << umis.size() << " barcodes extracted.")

    std::unordered_map<Umi, std::vector<size_t>> umi_to_reads;
    group_reads_by_umi(umis, umi_to_reads);

    report_umi_group_sizes(umi_to_reads);
    report_umi_length_distribution(umi_to_reads);

    INFO("Compressing barcodes.")
    std::vector<seqan::Dna5String> compressed_umis;
    compressed_umis.reserve(umis.size());
    for (const auto& umi : umis) {
        compressed_umis.push_back(compress_string(umi));
    }

    std::vector<seqan::CharString> output_ids;
    output_ids.reserve(input_ids.size());
    for (size_t i = 0; i < input_ids.size(); i ++) {
        std::string id = seqan_string_to_string(input_ids[i]);
        const std::string umi = seqan_string_to_string(umis[i]);
        const std::string compressed_umi = seqan_string_to_string(compressed_umis[i]);
        id.replace(id.find(umi), id.length(), compressed_umi);
        output_ids.push_back(seqan::CharString(id.c_str()));
    }
    write_output(params.output_file, output_ids, input_reads);

    std::unordered_map<Umi, std::vector<size_t>> compressed_umi_to_reads;
    group_reads_by_umi(compressed_umis, compressed_umi_to_reads);
    for (const auto& entry : compressed_umi_to_reads) {
        const auto ids = entry.second;
        if (ids.size() < 20) continue;
        const auto umi = entry.first.GetString();
        std::cout << ids.size() << "\n";
        for (size_t idx : ids) {
            std::cout << ">" << output_ids[idx] << "\n" << compress_string(input_reads[idx]) << "\n";
        }
        std::cout << std::endl;
    }

    report_umi_group_sizes(compressed_umi_to_reads);
    report_umi_length_distribution(compressed_umi_to_reads);

    return 0;
}
