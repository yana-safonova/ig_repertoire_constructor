#include <logger/logger.hpp>
#include <logger/log_writers.hpp>
#include <segfault_handler.hpp>
#include <perfcounter.hpp>
#include <boost/program_options.hpp>
#include <seqan/file.h>
#include <seqan/store.h>
#include <uchar.h>
#include <bitset>
#include "../ig_tools/utils/string_tools.hpp"

struct Options {
    struct PcrOptions {
        size_t cycles_count;
        double error_prob_first;
        double error_prob_last;
        double amplification_rate;
        double chimeras_rate;
        // 1 for barcode going with the left half, 2 for the right half, 3 for random choice
        size_t barcode_position;
    };

    std::string repertoire_file_path;
    std::string output_file_path;
    size_t barcode_length;
    PcrOptions pcr_options;
    size_t output_estimation_limit;
};

Options parse_options(int argc, const char* const* argv) {
    Options options;
    Options::PcrOptions& pcrOptions = options.pcr_options;
    namespace po = boost::program_options;
    po::options_description cmdline_options("Command line options");
    cmdline_options.add_options()
            ("input-file,i", po::value<std::string>(&options.repertoire_file_path)->required(),
             "name of the input repertoire file (FASTA|FASTQ)")
            ("output-file,o", po::value<std::string>(&options.output_file_path)->required(),
             "output file for final repertoire")
            ("umi-length,l", po::value<size_t>(&options.barcode_length)->default_value(14),
             "length of generated barcodes (defaults to 14)")
            ("pcr-cycles,c", po::value<size_t>(&pcrOptions.cycles_count)->default_value(25),
             "number of PCR cycles to simulate (defaults to 25)")
            ("pcr-error1,e", po::value<double>(&pcrOptions.error_prob_first)->required(),
             "probability of PCR error on the first PCR cycle")
            ("pcr-error2,E", po::value<double>(&pcrOptions.error_prob_last)->required(),
             "probability of PCR error on the last PCR cycle (probability on the other cycles are interpolated)")
            ("pcr-rate,r", po::value<double>(&pcrOptions.amplification_rate)->default_value(0.1),
             "probability for each molecule to be amplified on each PCR cycle (defaults to 0.1)")
            ("chimeras-rate,k", po::value<double>(&pcrOptions.chimeras_rate)->default_value(0.001),
             "percentage of chimeric reads appering on each cycle (defaults to 0.001)")
            ("output-limit,m", po::value<size_t>(&options.output_estimation_limit)->default_value(static_cast<size_t>(1e8)),
             "the program will exit if expected number of reads exceeds this parameter (defaults to 100,000,000)")
            ("barcode-position,b", po::value<size_t>(&pcrOptions.barcode_position)->default_value(3),
             "indicator of barcode position in the read, used for chimeras simulation. "
             "1 for barcode going with the left half, 2 for the right half, 3 for random choice (defaults to 3)")
            ;
    po::variables_map vm;
    store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
    po::notify(vm);
    return options;
}

void create_console_logger() {
    using namespace logging;
    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

std::vector<seqan::Dna5String> read_repertoire(Options options) {
    std::vector<seqan::CharString> read_ids;
    std::vector<seqan::Dna5String> reads;
    seqan::SeqFileIn reads_file(options.repertoire_file_path.c_str());
    readRecords(read_ids, reads, reads_file);
    double exp_reads_count = static_cast<double>(reads.size()) *
                             pow(1.0 + options.pcr_options.amplification_rate + options.pcr_options.chimeras_rate, static_cast<double>(options.pcr_options.cycles_count));
    VERIFY(exp_reads_count <= options.output_estimation_limit);
    return reads;
}

seqan::Dna5String generate_barcode(size_t length) {
    seqan::Dna5String barcode;
    for (size_t i = 0; i < length; i ++) {
        barcode += std::rand() % 4;
    }
    return barcode;
}

std::vector<seqan::Dna5String> generate_barcodes(size_t count, size_t barcode_length) {
    std::vector<seqan::Dna5String> barcodes(count);
    for (size_t i = 0; i < count; i ++) {
        barcodes[i] = generate_barcode(barcode_length);
    }
    return barcodes;
}

void amplify(std::vector<seqan::Dna5String>& reads, std::vector<seqan::Dna5String>& barcodes, std::vector<seqan::CharString>& ids, double pcr_error_prob,
             Options::PcrOptions options, std::minstd_rand0& random_engine) {
    size_t size = reads.size();
    for (size_t read_idx = 0; read_idx < size; read_idx ++) {
        if (std::rand() <= static_cast<double>(RAND_MAX) * options.amplification_rate) {
            seqan::Dna5String read = reads[read_idx];
            seqan::Dna5String barcode = barcodes[read_idx];
            for (size_t pos = 0; pos < length(barcode) + length(read); pos ++) {
                if (std::rand() <= static_cast<double>(RAND_MAX) * pcr_error_prob) {
                    auto current_value = pos < length(barcode) ? barcode[pos] : read[pos - length(barcode)];
                    int new_value_candidate = std::rand() % 3;
                    int new_value = new_value_candidate + (new_value_candidate >= current_value ? 1 : 0);
                    pos < length(barcode) ? barcode[pos] : read[pos - length(barcode)] = new_value;
                }
            }
            barcodes.push_back(barcode);
            reads.push_back(read);
            ids.emplace_back(std::to_string(reads.size()) + "_mutated_from_" + std::to_string(reads.size()));
        }
    }
    std::uniform_int_distribution<size_t> read_distribution(0, size - 1);
    std::uniform_int_distribution<size_t> barcode_position_distribution(1, std::bitset<2>(options.barcode_position).count());
    for (size_t chimera = 0; chimera < static_cast<double>(size) * options.chimeras_rate; chimera ++) {
        size_t left_idx = read_distribution(random_engine);
        size_t right_idx = read_distribution(random_engine);
        size_t barcode_idx = (options.barcode_position & 1) && barcode_position_distribution(random_engine) == 1 ? left_idx : right_idx;
        barcodes.push_back(barcodes[barcode_idx]);
        std::string left = seqan_string_to_string(reads[left_idx]);
        left = left.substr(0, left.length() / 2);
        std::string right = seqan_string_to_string(reads[right_idx]);
        right = right.substr(right.length() / 2);
        reads.push_back(seqan::Dna5String(left + right));
        ids.emplace_back(std::to_string(reads.size()) + "_chimera_from_" + std::to_string(left_idx) + "_" + std::to_string(right_idx));
    }
    VERIFY(reads.size() == barcodes.size() && reads.size() == ids.size());
}

void simulate_pcr(std::vector<seqan::Dna5String>& reads, std::vector<seqan::Dna5String>& barcodes, std::vector<seqan::CharString>& ids, Options::PcrOptions options, std::minstd_rand0& random_engine) {
    for (size_t i = 0; i < options.cycles_count; i ++) {
        double pcr_error_prob = options.error_prob_first + (options.error_prob_last - options.error_prob_first) * static_cast<double>(i) / static_cast<double>(options.cycles_count - 1);
        amplify(reads, barcodes, ids, pcr_error_prob, options, random_engine);
    }
}

void write_repertoire(const std::vector<seqan::Dna5String>& reads, const std::vector<seqan::Dna5String>& barcodes, const std::vector<seqan::CharString>& ids, std::vector<size_t> perm, std::string file_name) {
    seqan::SeqFileOut output_file(file_name.c_str());
    for (size_t i = 0; i < reads.size(); i ++) {
        std::stringstream sstr;
        sstr << ids[perm[i]] << "_UMI:" << barcodes[perm[i]];
        seqan::writeRecord(output_file, sstr.str(), reads[perm[i]]);
    }
}

int main(int argc, const char* const* argv) {
    perf_counter pc;
    segfault_handler sh;
    create_console_logger();
    std::srand(495634);
    auto random_engine = std::default_random_engine(8356);

    const Options& options = parse_options(argc, argv);

    auto reads = read_repertoire(options);
    std::vector<seqan::CharString> ids(reads.size());
    for (size_t i = 0; i < reads.size(); i ++) {
        ids[i] = "original_" + std::to_string(i);
    }
    auto barcodes = generate_barcodes(reads.size(), options.barcode_length);

    simulate_pcr(reads, barcodes, ids, options.pcr_options, random_engine);

    std::vector<size_t> perm(reads.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), random_engine);

    write_repertoire(reads, barcodes, ids, perm, options.output_file_path);
}
