#include <string>
#include <logger/log_writers.hpp>
#include <segfault_handler.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include <ion_utils.hpp>
#include "../../umi_experiments/utils.hpp"

namespace fs = boost::filesystem;

namespace {
    struct Params {
        std::string input_germline_dir;
        std::string output_germline_dir;
    };

    bool read_args(int argc, const char *const *argv, Params& params) {
        namespace po = boost::program_options;
        po::options_description cmdl_options;
        cmdl_options.add_options()
                ("help,h", "print help message")
                ("input,i", po::value<std::string>(&params.input_germline_dir)->required(), "input file with reads")
                ("output,o", po::value<std::string>(&params.output_germline_dir)->required(), "output file with compressed barcodes");
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(cmdl_options).run(), vm);
        if (vm.count("help") || argc == 1) {
            std::cout << cmdl_options << std::endl;
            return false;
        }
        po::notify(vm);
        return true;
    }

    void compress_file(const fs::path& from, const fs::path& to) {
        seqan::SeqFileIn reads_file(from.c_str());
        std::vector<seqan::CharString> ids;
        std::vector<seqan::Dna5String> seqs;
        readRecords(ids, seqs, reads_file);
        std::vector<seqan::Dna5String> compressed;
        compressed.reserve(seqs.size());
        for (const auto& seq : seqs) {
            compressed.push_back(compress_string(seq));
        }
        seqan::SeqFileOut output_file(to.c_str());
        writeRecords(output_file, ids, compressed);
    }

    void compress_recursively(const fs::path& from_dir, const fs::path& to_dir) {
        INFO("Compressing files from " << from_dir << " to " << to_dir);
        for (fs::directory_iterator itr(from_dir); itr != fs::directory_iterator(); ++ itr) {
            const auto& path = itr->path();
            auto to_child = to_dir;
            to_child /= path.filename();
            if (fs::is_directory(path)) {
                fs::create_directory(to_child);
                compress_recursively(path, to_child);
            } else if (fs::is_regular_file(path)) {
                INFO("Compressing " << path << " to " << to_child);
                compress_file(path, to_child);
            } else {
                ERROR("Unknown directory entry " << path);
            }
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

    INFO("Creating compressed germline from " << params.input_germline_dir << " at " << params.output_germline_dir);
    if (fs::exists(params.output_germline_dir)) {
        INFO("Removing existing output directory " << params.output_germline_dir);
        if (!fs::remove_all(params.output_germline_dir)) {
            FATAL_ERROR("Could not remove existing output directory");
        }
    }
    fs::create_directory(params.output_germline_dir);

    compress_recursively(params.input_germline_dir, params.output_germline_dir);

    return 0;
}
