#include <logger/logger.hpp>
#include <logger/log_writers.hpp>
#include <verify.hpp>
#include <segfault_handler.hpp>

#include <read_archive.hpp>
#include <copy_file.hpp>
#include <boost/program_options.hpp>

#include <seqan/seq_io.h>

namespace po = boost::program_options;

// void prepare_output_dir(const vj_finder::VJFinderConfig::IOParams::OutputParams::OutputFiles & of) {
//     path::make_dir(of.output_dir);
// }

po::variables_map parse_command_line(int argc, char** argv) {
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("input,i", po::value<std::string>(), "name of the input FASTQ file")
        ("output,o", po::value<std::string>(), "name of the output FASTQ file")
        ("threshold,t", po::value<char>()->default_value(char(2)),
         "quality threshold: trim tails with quality <= threshold")
        ;

    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc, (const char **) argv).options(desc).run();
    po::store(parsed, vm);
    vm.notify();

    if (vm.count("help")) {
        std::cout << desc;
        exit(0);
    }

    return vm;
}

int main(int argc, char **argv) {
    omp_set_num_threads(1);

    segfault_handler sh;
    perf_counter pc;
    po::variables_map vm = parse_command_line(argc, argv);
    std::string input = vm["input"].as<std::string>();
    std::string output = vm["output"].as<std::string>();
    char threshold = vm["threshold"].as<char>() + '!';

    seqan::SeqFileIn seqan_input(seqan::toCString(input));
    seqan::SeqFileOut seqan_output(seqan::toCString(output));
    // seqan::CharString id;
    // seqan::Dna5String seq;
    std::string id, seq, qual;

    size_t i = 0;
    while(!seqan::atEnd(seqan_input)) {
        if (not (i % 10000))
            std::cout << "Processed " << i << " sequences\n";

        seqan::readRecord(id, seq, qual, seqan_input);
        size_t left = 0, right = qual.length();
        while (left < qual.length() and qual[left] <= threshold)
            left++;

        if (left == qual.length()) {
            seqan::writeRecord(seqan_output, id, seq, qual);
            continue;
        }

        while (right and qual[right - 1] <= threshold)
            right--;

        seqan::writeRecord(seqan_output, id,
                           seq.substr(left, right - left),
                           qual.substr(left, right - left));
        i++;
    }

    return 0;
}
