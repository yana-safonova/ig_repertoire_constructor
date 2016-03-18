#include <boost/filesystem.hpp>
#include <segfault_handler.hpp>
#include <logger/log_writers.hpp>
#include "umi_utils.hpp"
#include "utils.hpp"

bool parse_cmdline(int argc, char **argv, std::string& input_file, std::string& output_dir) {
    if (argc != 3) {
        std::cout << "Usage: <input file> <output dir>" << std::endl;
        return false;
    }
    input_file = argv[1];
    output_dir = argv[2];
    return true;
}

void print_reads_by_umi(std::string& output_dir_name, std::vector<seqan::CharString>& input_ids, std::vector<seqan::Dna5String>& input_reads,
                        std::vector<seqan::Dna5String>& input_umi, bool group_by_size, size_t group_size_threshold) {
    INFO("Grouping reads by UMI");
    std::unordered_map<Umi, std::vector<size_t>> umi_to_reads;
    for (size_t i = 0; i < input_ids.size(); i ++) {
        umi_to_reads[Umi(input_umi[i])].push_back(i);
    }

    auto out_dir = boost::filesystem::absolute(output_dir_name);
    INFO("Removing " << out_dir << " and all its contents");
    boost::filesystem::remove_all(out_dir);
    boost::filesystem::create_directory(out_dir);

    INFO("Storing result in " << out_dir);
    for (auto entry : umi_to_reads) {
        auto umi = entry.first;
        auto reads = entry.second;
        size_t reads_count = reads.size();
        if (reads_count < group_size_threshold) continue;
        auto size_dir = group_by_size ? boost::filesystem::path(out_dir).append(std::to_string(reads_count)) : out_dir;
        if (!boost::filesystem::exists(size_dir)) {
            boost::filesystem::create_directory(size_dir);
        }
        const auto umi_file_path = boost::filesystem::path(size_dir).append(seqan_string_to_string(umi.GetString())).replace_extension(".fasta");
        seqan::SeqFileOut output_file(umi_file_path.c_str());
        std::vector<seqan::CharString> umi_read_ids(reads_count);
        std::vector<seqan::Dna5String> umi_reads(reads_count);
        for (size_t i = 0; i < reads_count; i ++) {
            umi_read_ids[i] = input_ids[reads[i]];
            umi_reads[i] = input_reads[reads[i]];
        }
        writeRecords(output_file, umi_read_ids, umi_reads);
    }
    INFO("Result written to " << out_dir);
}

int main(int argc, char** argv) {
    segfault_handler sh;
    create_console_logger();

    std::string input_file;
    std::string output_dir;
    if (!parse_cmdline(argc, argv, input_file, output_dir)) return 0;

    INFO("Reading fastq");
    seqan::SeqFileIn seqFileIn_input(input_file.c_str());
    std::vector<seqan::CharString> input_ids;
    std::vector<seqan::Dna5String> input_reads;
    std::vector<seqan::DnaQString> input_qual;
    readRecords(input_ids, input_reads, input_qual, seqFileIn_input);

    INFO("Extracting UMI data");
    std::vector<seqan::Dna5String> input_umi;
    std::vector<seqan::DnaQString> input_umi_qual;
    extract_barcodes_from_read_ids(input_ids, input_umi, input_umi_qual);

    INFO("Printing reads by UMI");
    print_reads_by_umi(output_dir, input_ids, input_reads, input_umi, true, 5);
}
