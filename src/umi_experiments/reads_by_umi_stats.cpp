#include <atomic>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <segfault_handler.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include "utils.hpp"
#include "umi_utils.hpp"
#include "dist_distribution_stats.hpp"

bool read_args(int argc, char **argv, std::string& reads_file, std::string& output_dir, unsigned& thread_count) {
    namespace po = boost::program_options;
    po::options_description cmdl_options("Is this needed?");
    cmdl_options.add_options()
            ("help,h", "print help message")
            ("reads,r", po::value<std::string>(&reads_file)->required(), "input file with reads")
            ("output,o", po::value<std::string>(&output_dir)->default_value(""), "output directory")
            ("threads,t", po::value<unsigned>(&thread_count)->default_value(1), "number of running threads")
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

void print_umi_reads_distribution_by_size(unordered_map<Umi, vector<size_t>>& umi_to_reads) {
    std::map<size_t, size_t> size_to_count;
    for (auto& entry : umi_to_reads) {
        size_to_count[entry.second.size()] ++;
    }
    for (auto& entry : size_to_count) {
        std::cout << "size " << entry.first << ", count " << entry.second << std::endl;
    }
}

int main(int argc, char **argv) {
    segfault_handler sh;
    create_console_logger();
    std::string reads_file;
    std::string output_dir;
    unsigned thread_count;
    if (!read_args(argc, argv, reads_file, output_dir, thread_count)) {
        return 0;
    }

    INFO("Reading fastq");
    seqan::SeqFileIn seqFileIn_input(reads_file.c_str());
    std::vector<seqan::CharString> input_ids;
    std::vector<seqan::Dna5String> input_reads;
    readRecords(input_ids, input_reads, seqFileIn_input);
    INFO(input_ids.size() << " records read");

    INFO("Extracting barcodes");
    std::vector<seqan::Dna5String> umis;
    {
        std::vector<seqan::DnaQString> umi_quals;
        extract_barcodes_from_read_ids(input_ids, umis, umi_quals);
    }

    std::unordered_map<Umi, std::vector<size_t> > umi_to_reads;
    group_reads_by_umi(umis, umi_to_reads);
    INFO("stats by sizes:");
    print_umi_reads_distribution_by_size(umi_to_reads);

    INFO("Calculating stats");
    auto stats = DistDistributionStats::GetStats(input_reads, umi_to_reads, thread_count);

    if (output_dir != "") {
        INFO("Saving results");
        if (!boost::filesystem::exists(output_dir)) {
            INFO("Creating directory " << output_dir);
            boost::filesystem::create_directory(output_dir);
        }
        for (auto umi_size : stats.GetSizes()) {
            const auto file_path = boost::filesystem::path(output_dir).append(std::to_string(umi_size));
            std::ofstream output(file_path.string());
            output << stats.ToString(umi_size) << std::endl;
        }
        INFO("Stats printed to " << output_dir);
    } else {
        for (auto umi_size : stats.GetSizes()) {
            std::cout << stats.ToString(umi_size);
        }
    }
}
