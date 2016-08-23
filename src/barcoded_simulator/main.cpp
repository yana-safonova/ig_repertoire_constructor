#include <logger/logger.hpp>
#include <logger/log_writers.hpp>
#include <segfault_handler.hpp>
#include <perfcounter.hpp>
#include <boost/program_options.hpp>
#include <seqan/file.h>
#include <seqan/store.h>

struct Options {
    std::string repertoire_file_path;
    std::string output_file_path;
    size_t barcode_length;
    size_t pcr_cycles;
    float pcr_error_prob_first;
    float pcr_error_prob_last;
    float amplification_rate;
};

Options parse_options(int argc, const char* const* argv) {
    Options options;
    namespace po = boost::program_options;
    po::options_description cmdline_options("Command line options");
    cmdline_options.add_options()
            ("input-file,i", po::value<std::string>(&options.repertoire_file_path)->required(),
             "name of the input repertoire file (FASTA|FASTQ)")
            ("output-file,o", po::value<std::string>(&options.output_file_path)->required(),
             "output file for final repertoire")
            ("umi-length,l", po::value<size_t>(&options.barcode_length)->default_value(14),
             "length of generated barcodes (defaults to 14)")
            ("pcr-cycles,c", po::value<size_t>(&options.pcr_cycles)->default_value(25),
             "number of PCR cycles to simulate (defaults to 25)")
            ("pcr-error1,e", po::value<float>(&options.pcr_error_prob_first),
             "probability of PCR error on the first PCR cycle")
            ("pcr-error2,E", po::value<float>(&options.pcr_error_prob_last),
             "probability of PCR error on the last PCR cycle (probability on the other cycles are interpolated)")
            ("pcr-rate,r", po::value<float>(&options.amplification_rate)->default_value(1.0),
             "probability for each molecule to be amplified on each PCR cycle (defaults to 1.0)")
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

std::vector<seqan::Dna5String> read_repertoire(std::string file_name) {
    std::vector<seqan::CharString> read_ids;
    std::vector<seqan::Dna5String> reads;
    seqan::SeqFileIn reads_file(file_name.c_str());
    readRecords(read_ids, reads, reads_file);
    return reads;
}

seqan::Dna5String generate_barcode(size_t length) {
    seqan::Dna5String barcode;
    for (size_t i = 0; i < length; i ++) {
        barcode += std::rand() % 4;
    }
    return barcode;
}

std::vector<seqan::Dna5String> attach_barcodes(size_t count, size_t barcode_length) {
    std::vector<seqan::Dna5String> barcodes(count);
    for (size_t i = 0; i < count; i ++) {
        barcodes[i] = generate_barcode(barcode_length);
    }
    return barcodes;
}

void amplify(std::vector<seqan::Dna5String> reads, std::vector<seqan::Dna5String> barcodes, float pcr_error_prob, float amplification_rate) {
    size_t size = reads.size();
    for (size_t read_idx = 0; read_idx < size; read_idx ++) {
        if (std::rand() <= static_cast<float>(RAND_MAX) * amplification_rate) {
            seqan::Dna5String read = reads[read_idx];
            seqan::Dna5String barcode = barcodes[read_idx];
            for (size_t pos = 0; pos < length(barcode) + length(read); pos ++) {
                if (std::rand() <= RAND_MAX * pcr_error_prob) {
                    auto current_value = pos < length(barcode) ? barcode[pos] : read[pos - length(barcode)];
                    int new_value_candidate = std::rand() % 3;
                    int new_value = new_value_candidate + (new_value_candidate >= current_value ? 1 : 0);
                    pos < length(barcode) ? barcode[pos] : read[pos - length(barcode)] = new_value;
                }
            }
            barcodes.push_back(barcode);
            reads.push_back(read);
        }
    }
}

void simulate_pcr(std::vector<seqan::Dna5String>& reads, std::vector<seqan::Dna5String>& barcodes, size_t cycles,
                  float pcr_error_prob_first, float pcr_error_prob_last, float amplification_rate) {
    for (size_t i = 0; i < cycles; i ++) {
        amplify(reads, barcodes, pcr_error_prob_first + (pcr_error_prob_last - pcr_error_prob_first) * static_cast<float>(i) / (cycles - 1), amplification_rate);
    }
}

void write_repertoire(std::vector<seqan::Dna5String> reads, std::vector<seqan::Dna5String> barcodes, std::vector<size_t> perm, std::string file_name) {
    seqan::SeqFileOut output_file(file_name.c_str());
    std::vector<seqan::CharString> ids;
    for (size_t i = 0; i < reads.size(); i ++) {
        std::stringstream sstr;
        sstr << "read_" << i << "_UMI:" << barcodes[perm[i]];
        seqan::writeRecord(output_file, sstr.str(), reads[perm[i]]);
    }
}

int main(int argc, const char* const* argv) {
    perf_counter pc;
    segfault_handler sh;
    create_console_logger();
    const Options& options = parse_options(argc, argv);

    auto reads = read_repertoire(options.repertoire_file_path);
    std::srand(495634);
    auto barcodes = attach_barcodes(reads.size(), options.barcode_length);
    simulate_pcr(reads, barcodes, options.pcr_cycles, options.pcr_error_prob_first, options.pcr_error_prob_last, options.amplification_rate);
    std::vector<size_t> perm(reads.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), std::default_random_engine(8972356));
    write_repertoire(reads, barcodes, perm, options.output_file_path);
}
