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
    const size_t MIN_READ_LENGTH = 50;

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
}


int main(int argc, const char* const* argv) {
    segfault_handler sh;
    create_console_logger(logging::L_TRACE);

    Params params;
    if (!read_args(argc, argv, params)) {
        return 1;
    }

    std::vector<seqan::CharString> input_ids;
    std::vector<seqan::Dna5String> input_reads;
    read_input(params.input_file, input_ids, input_reads);

    std::vector<seqan::Dna5String> umis;
    extract_barcodes_from_read_ids(input_ids, umis);
    INFO("Total " << umis.size() << " barcodes extracted.")

    INFO("Compressing barcodes and reads.")
    std::vector<seqan::Dna5String> compressed_umis;
    std::vector<seqan::Dna5String> compressed_reads;
    std::vector<seqan::CharString> output_ids;
    compressed_umis.reserve(umis.size());
    compressed_reads.reserve(umis.size());
    output_ids.reserve(input_ids.size());
    for (size_t i = 0; i < umis.size(); i ++) {
        const auto read = compress_string(input_reads[i]);
        if (length(read) >= MIN_READ_LENGTH) {
            const auto umi = compress_string(umis[i]);
            compressed_umis.push_back(umi);
            compressed_reads.push_back(read);

            std::string id = seqan_string_to_string(input_ids[i]);
            const std::string umi_str = seqan_string_to_string(umis[i]);
            const std::string compressed_umi_str = seqan_string_to_string(umi);
            id.replace(id.length() - umi_str.length(), umi_str.length(), compressed_umi_str);
            output_ids.emplace_back(id);
        }
    }
    INFO((umis.size() - compressed_umis.size()) << " reads filtered out as too short (shorter than " << MIN_READ_LENGTH << ").");

    write_output(params.output_file, output_ids, compressed_reads);

    return 0;
}
