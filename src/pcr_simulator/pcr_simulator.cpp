#include <boost/algorithm/string.hpp>
#include <verify.hpp>
#include <logger/logger.hpp>
#include <bitset>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include "pcr_simulator.hpp"
#include "../ig_tools/utils/string_tools.hpp"

std::default_random_engine PcrSimulator::random_engine_;

void PcrSimulator::read_repertoire(const std::string& repertoire_file_path) {
    seqan::SeqFileIn reads_file(repertoire_file_path.c_str());
    readRecords(original_ids_, original_reads_, reads_file);
}

void PcrSimulator::amplify(const size_t output_estimation_limit) {
    check_limit(output_estimation_limit);

    size_t current = 0;
    for (size_t i = 0; i < original_ids_.size(); i ++) {
        auto id = seqan_string_to_string(original_ids_[i]);
        if (boost::algorithm::ends_with(id, "_copy_1")) {
            current ++;
        }
        read_to_compressed_.push_back(current - 1);
    }

    amplified_reads_ = std::vector<seqan::Dna5String>(original_reads_);
    amplified_ids_ = std::vector<seqan::CharString>(original_reads_.size());
    for (size_t i = 0; i < amplified_ids_.size(); i ++) {
        amplified_ids_[i] = "original_" + std::to_string(i);
    }
    generate_barcodes();

    simulate_pcr();

    perm_ = std::vector<size_t>(amplified_reads_.size());
    std::iota(perm_.begin(), perm_.end(), 0);
    std::shuffle(perm_.begin(), perm_.end(), random_engine_);
}

void PcrSimulator::check_limit(const size_t output_estimation_limit) {
    double exp_reads_count = static_cast<double>(original_reads_.size()) *
                             pow(1.0 + options_.amplification_rate + options_.chimeras_rate, static_cast<double>(options_.cycles_count));
    VERIFY(exp_reads_count <= output_estimation_limit);
}

seqan::Dna5String PcrSimulator::generate_barcode() {
    seqan::Dna5String barcode;
    for (size_t i = 0; i < options_.barcode_length; i ++) {
        barcode += std::rand() % 4;
    }
    return barcode;
}

void PcrSimulator::generate_barcodes() {
    amplified_barcodes_ = std::vector<seqan::Dna5String>(amplified_reads_.size());
    for (auto& barcode : amplified_barcodes_) {
        barcode = generate_barcode();
    }
}

void PcrSimulator::report_average_error_rate() {
    size_t total_errors = 0;
    size_t interesting_reads = 0;
    for (size_t errors : read_error_count_) {
        if (errors < std::numeric_limits<size_t>::max()) {
            total_errors += errors;
            interesting_reads ++;
        }
    }
    INFO("Average amount of errors per read is " << ((double) total_errors / (double) interesting_reads) << " (total " << total_errors << " in " << interesting_reads << " of " << read_error_count_.size() << " reads)");
    INFO("Average amount of errors per barcode is " << ((double) barcode_error_count_ / (double) read_error_count_.size()) << " (total " << barcode_error_count_ << " in " << read_error_count_.size() << " reads)");
}

void PcrSimulator::amplify_sequences(double pcr_error_prob) {
    size_t size = amplified_reads_.size();
    for (size_t read_idx = 0; read_idx < size; read_idx ++) {
        if (std::rand() <= static_cast<double>(RAND_MAX) * options_.amplification_rate) {
            seqan::Dna5String read = amplified_reads_[read_idx];
            seqan::Dna5String barcode = amplified_barcodes_[read_idx];
            size_t errors = read_error_count_[read_idx];
            for (size_t pos = 0; pos < length(barcode) + length(read); pos ++) {
                if (std::rand() <= static_cast<double>(RAND_MAX) * pcr_error_prob) {
                    auto current_value = pos < length(barcode) ? barcode[pos] : read[pos - length(barcode)];
                    int new_value_candidate = std::rand() % 3;
                    int new_value = new_value_candidate + (new_value_candidate >= current_value ? 1 : 0);
                    pos < length(barcode) ? barcode[pos] : read[pos - length(barcode)] = new_value;
                    if (pos >= length(barcode)) {
                        errors ++;
                    } else {
                        barcode_error_count_ ++;
                    }
                }
            }
            amplified_barcodes_.push_back(barcode);
            amplified_reads_.push_back(read);
            amplified_ids_.emplace_back(std::to_string(amplified_reads_.size()) +
                             (seqan_string_to_string(amplified_ids_[read_idx]).find("chimera") == std::string::npos ? "" : "_chimera") +
                             "_mutated_from_" +
                             std::to_string(read_idx));
            read_error_count_.push_back(errors);
            read_to_compressed_.push_back(read_to_compressed_[read_idx]);
        }
    }
    std::uniform_int_distribution<size_t> read_distribution(0, size - 1);
    std::uniform_int_distribution<size_t> barcode_position_distribution(1, std::bitset<2>(options_.barcode_position).count());
    for (size_t chimera = 0; chimera < static_cast<double>(size) * options_.chimeras_rate; chimera ++) {
        size_t left_idx = read_distribution(random_engine_);
        size_t right_idx = read_distribution(random_engine_);
        size_t barcode_idx = (options_.barcode_position & 1) && barcode_position_distribution(random_engine_) == 1 ? left_idx : right_idx;
        amplified_barcodes_.push_back(amplified_barcodes_[barcode_idx]);
        std::string left = seqan_string_to_string(amplified_reads_[left_idx]);
        left = left.substr(0, left.length() / 2);
        std::string right = seqan_string_to_string(amplified_reads_[right_idx]);
        right = right.substr(right.length() / 2);
        amplified_reads_.push_back(seqan::Dna5String(left + right));
        read_error_count_.push_back(std::numeric_limits<size_t>::max());
        read_to_compressed_.push_back(std::numeric_limits<size_t>::max());
        amplified_ids_.emplace_back(std::to_string(amplified_reads_.size()) + "_chimera_from_" + std::to_string(left_idx) + "_" + std::to_string(right_idx));
    }
    VERIFY(amplified_reads_.size() == amplified_barcodes_.size() &&
                   amplified_reads_.size() == amplified_ids_.size() &&
                   amplified_reads_.size() == read_to_compressed_.size() &&
                   amplified_reads_.size() == read_error_count_.size());

//    report_average_error_rate(read_error_count);
}

void PcrSimulator::simulate_pcr() {
    read_error_count_ = std::vector<size_t>(amplified_reads_.size(), 0);
    barcode_error_count_ = 0;
    for (size_t i = 0; i < options_.cycles_count; i ++) {
        double pcr_error_prob = options_.error_prob_first + (options_.error_prob_last - options_.error_prob_first) * static_cast<double>(i) / static_cast<double>(options_.cycles_count - 1);
        amplify_sequences(pcr_error_prob);
    }

    report_average_error_rate();
}

void PcrSimulator::write_results(const std::string& output_dir_path) {
    boost::filesystem::create_directory(output_dir_path);
    write_repertoire(boost::filesystem::path(output_dir_path).append("amplified.fasta").string());
    write_compressed(boost::filesystem::path(output_dir_path).append("repertoire_comp.fasta").string());
}

void PcrSimulator::write_repertoire(const std::string& path) {
    seqan::SeqFileOut output_file(path.c_str());
    for (size_t i = 0; i < amplified_reads_.size(); i ++) {
        std::stringstream sstr;
        sstr << amplified_ids_[perm_[i]] << "_UMI:" << amplified_barcodes_[perm_[i]];
        seqan::writeRecord(output_file, sstr.str(), amplified_reads_[perm_[i]]);
    }
}

void PcrSimulator::write_compressed(const std::string& path) {
    std::map<size_t, size_t> id_to_count;
    for (size_t id : read_to_compressed_) {
        id_to_count[id] ++;
    }

    std::vector<seqan::CharString> compressed_ids;
    std::vector<seqan::Dna5String> compressed_reads;
    size_t current = 0;
    for (size_t i = 0; i < amplified_ids_.size(); i ++) {
        read_to_compressed_.push_back(current);
        auto id = seqan_string_to_string(amplified_ids_[i]);
        if (boost::algorithm::ends_with(id, "_copy_1")) {
            current ++;
            compressed_ids.emplace_back((boost::format("cluster___%d___size___%d") % (current - 1) % id_to_count[current - 1]).str());
            compressed_reads.push_back(amplified_reads_[i]);
        }
    }
    seqan::SeqFileOut compressed_file(path.c_str());
    seqan::writeRecords(compressed_file, compressed_ids, compressed_reads);
}
