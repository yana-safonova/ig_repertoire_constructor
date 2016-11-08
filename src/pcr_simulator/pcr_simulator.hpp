#pragma once

#include <seqan/file.h>
#include "options.hpp"

class PcrSimulator {
private:
    static std::default_random_engine random_engine_;

    const Options::SimulationOptions options_;
    std::vector<seqan::CharString> original_ids_;
    std::vector<seqan::Dna5String> original_reads_;
    std::vector<seqan::CharString> amplified_ids_;
    std::vector<seqan::Dna5String> amplified_reads_;
    std::vector<seqan::Dna5String> amplified_barcodes_;
    std::vector<size_t> perm_;
    std::vector<size_t> read_to_compressed_;
    std::vector<size_t> read_error_count_;
    size_t barcode_error_count_;

public:
    PcrSimulator(const Options::SimulationOptions& options) : options_(options) {
        random_engine_.seed(8356);
    };
    void read_repertoire(const std::string& repertoire_file_path);
    void amplify(size_t output_estimation_limit);
    void write_results(const std::string& output_dir_path);

private:
    void check_limit(size_t output_estimation_limit);
    seqan::Dna5String generate_barcode();
    void generate_barcodes();
    void report_average_error_rate();
    void simulate_pcr();
    void amplify_sequences(double pcr_error_prob);
    void write_repertoire(const std::string& path);
    void write_compressed(const std::string& path);
};
