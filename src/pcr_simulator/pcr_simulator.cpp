#include <boost/algorithm/string.hpp>
#include <verify.hpp>
#include <logger/logger.hpp>
#include <bitset>
#include <unordered_set>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include "pcr_simulator.hpp"
#include "../ig_tools/utils/string_tools.hpp"
#include "../umi_experiments/umi_utils.hpp"

std::default_random_engine PcrSimulator::random_engine_;
const size_t PcrSimulator::CHIMERIC_READ_ERROR_COUNT = std::numeric_limits<size_t>::max() / 2;
const size_t PcrSimulator::MAP_NOWHERE = std::numeric_limits<size_t>::max();

void PcrSimulator::ReadRepertoire(const std::string& repertoire_file_path) {
    seqan::SeqFileIn reads_file(repertoire_file_path.c_str());
    readRecords(original_ids_, original_reads_, reads_file);
    INFO("Total " << original_ids_.size() << " records read.");
    size_t total_length = 0;
    for (const auto& read : original_reads_) {
        total_length += length(read);
    }
    INFO("Average read length: " << static_cast<double>(total_length) / static_cast<double>(original_reads_.size()));
}

void PcrSimulator::Amplify(size_t output_estimation_limit) {
    CheckLimit(output_estimation_limit);

    size_t current = 0;
    for (size_t i = 0; i < original_ids_.size(); i ++) {
        auto id = seqan_string_to_string(original_ids_[i]);
        if (boost::algorithm::ends_with(id, "_copy_1")) {
            current ++;
        }
        read_to_compressed_.push_back(current - 1);
        read_to_original_.push_back(i);
    }

    amplified_reads_ = std::vector<Record>();
    for (size_t i = 0; i < original_reads_.size(); i ++) {
        amplified_reads_.emplace_back("original_" + std::to_string(i), original_reads_[i], GenerateBarcode(), 0);
    }

    {
        std::unordered_set<seqan::Dna5String> barcodes(amplified_reads_.size());
        for (const auto& read : amplified_reads_) {
            barcodes.insert(read.barcode);
        }
        INFO("Total " << barcodes.size() << " unique barcodes generated.");
    }

    SimulatePcr();

    perm_ = std::vector<size_t>(amplified_reads_.size());
    std::iota(perm_.begin(), perm_.end(), 0);
    std::shuffle(perm_.begin(), perm_.end(), random_engine_);
}

void PcrSimulator::CheckLimit(size_t output_estimation_limit) {
    double exp_reads_count = static_cast<double>(original_reads_.size()) *
                             pow(1.0 + options_.amplification_rate + options_.chimeras_rate, static_cast<double>(options_.cycles_count));
    VERIFY(exp_reads_count <= (double) output_estimation_limit);
}

seqan::Dna5String PcrSimulator::GenerateBarcode() {
    seqan::Dna5String barcode;
    for (size_t i = 0; i < options_.barcode_length; i ++) {
        barcode += std::rand() % 4;
    }
    return barcode;
}

void PcrSimulator::ReportAverageErrorRate() {
    size_t total_errors = 0;
    size_t interesting_reads = 0;
    for (const auto& record : amplified_reads_) {
        if (record.error_count < CHIMERIC_READ_ERROR_COUNT) {
            total_errors += record.error_count;
            interesting_reads ++;
        }
    }
    INFO("Average amount of errors per read is " << ((double) total_errors / (double) interesting_reads) << " (total " << total_errors << " in " << interesting_reads << " of " << amplified_reads_.size() << " reads)");
    INFO("Average amount of errors per barcode is " << ((double) barcode_error_count_ / (double) amplified_reads_.size()) << " (total " << barcode_error_count_ << " in " << amplified_reads_.size() << " reads)");
}

void PcrSimulator::AmplifySequences(double pcr_error_prob) {
    size_t size = amplified_reads_.size();
    for (size_t read_idx = 0; read_idx < size; read_idx ++) {
        if (std::rand() <= static_cast<double>(RAND_MAX) * options_.amplification_rate) {
            seqan::Dna5String new_read = amplified_reads_[read_idx].read;
            seqan::Dna5String new_barcode = amplified_reads_[read_idx].barcode;
            size_t errors = amplified_reads_[read_idx].error_count;
            for (size_t pos = 0; pos < length(new_barcode) + length(new_read); pos ++) {
                if (std::rand() <= static_cast<double>(RAND_MAX) * pcr_error_prob) {
                    auto current_value = pos < length(new_barcode) ? new_barcode[pos] : new_read[pos - length(new_barcode)];
                    int new_value_candidate = std::rand() % 3;
                    int new_value = new_value_candidate + (new_value_candidate >= current_value ? 1 : 0);
                    if (pos >= length(new_barcode)) {
                        new_read[pos - length(new_barcode)] = new_value;
                        errors ++;
                    } else {
                        new_barcode[pos] = new_value;
                        barcode_error_count_ ++;
                    }
                }
            }
            const std::string new_id = std::string(seqan_string_to_string(amplified_reads_[read_idx].id).find("chimera") == std::string::npos ? "" : "_chimera") +
                                       "_mutated_from_" + std::to_string(read_idx);
            AddRecord(new_id, new_read, new_barcode, errors);
            read_to_compressed_.push_back(read_to_compressed_[read_idx]);
            read_to_original_.push_back(read_to_original_[read_idx]);
        }
    }
    std::uniform_int_distribution<size_t> read_distribution(0, size - 1);
    std::uniform_int_distribution<size_t> barcode_position_distribution(1, std::bitset<2>(options_.barcode_position).count());
    for (size_t chimera = 0; chimera < (size_t) ((double) size * options_.chimeras_rate); chimera ++) {
        size_t left_idx = read_distribution(random_engine_);
        size_t right_idx = read_distribution(random_engine_);
        std::string left = seqan_string_to_string(amplified_reads_[left_idx].read);
        left = left.substr(0, left.length() / 2);
        std::string right = seqan_string_to_string(amplified_reads_[right_idx].read);
        right = right.substr(right.length() / 2);
        const std::string new_id = "_chimera_from_" + std::to_string(left_idx) + "_" + std::to_string(right_idx);
        const seqan::Dna5String new_read(left + right);
        size_t barcode_idx = (options_.barcode_position & 1) && barcode_position_distribution(random_engine_) == 1 ? left_idx : right_idx;
        const seqan::Dna5String& new_barcode = amplified_reads_[barcode_idx].barcode;
        AddChimericRecord(new_id, new_read, new_barcode);
        read_to_compressed_.push_back(MAP_NOWHERE);
        read_to_original_.push_back(MAP_NOWHERE);
    }
    VERIFY(amplified_reads_.size() == read_to_compressed_.size() && amplified_reads_.size() == read_to_original_.size());

//    ReportAverageErrorRate(read_error_count);
}

void PcrSimulator::AddRecord(const string& id, const seqan::Dna5String& read, const seqan::Dna5String& barcode, size_t error_count) {
    amplified_reads_.emplace_back(std::to_string(amplified_reads_.size()) + id, read, barcode, error_count);
}

void PcrSimulator::AddChimericRecord(const std::string& id, const seqan::Dna5String& read, const seqan::Dna5String& barcode){
    AddRecord(id, read, barcode, CHIMERIC_READ_ERROR_COUNT);
}

void PcrSimulator::SimulatePcr() {
    barcode_error_count_ = 0;
    double pcr_error_prob = options_.error_prob_first;
    for (size_t i = 0; i < options_.cycles_count; i ++) {
        AmplifySequences(pcr_error_prob);
        pcr_error_prob += (options_.error_prob_last - options_.error_prob_first) / static_cast<double>(options_.cycles_count - 1);
    }

    ReportAverageErrorRate();
}

void PcrSimulator::WriteResults(const std::string& output_dir_path) {
    boost::filesystem::create_directory(output_dir_path);
    UpdateReadIds();
    WriteRepertoire(boost::filesystem::path(output_dir_path).append("amplified.fasta").string());
    WriteRcm(boost::filesystem::path(output_dir_path).append("amplified_to_comp.rcm").string(), amplified_reads_, read_to_compressed_);
    WriteRcm(boost::filesystem::path(output_dir_path).append("amplified_to_orig.rcm").string(), amplified_reads_, read_to_original_);
    WriteCompressed(boost::filesystem::path(output_dir_path).append("repertoire_comp.fasta").string());
}

void PcrSimulator::UpdateReadIds() {
    for (auto& read : amplified_reads_) {
        std::stringstream ss;
        ss << read.id << "_UMI:" << read.barcode;
        read.id = seqan::CharString(ss.str());
    }
}

void PcrSimulator::WriteRepertoire(const std::string& path) {
    seqan::SeqFileOut output_file(path.c_str());
    for (size_t i = 0; i < amplified_reads_.size(); i ++) {
        std::stringstream sstr;
        seqan::writeRecord(output_file, amplified_reads_[perm_[i]].id, amplified_reads_[perm_[i]].read);
    }
}

void PcrSimulator::WriteCompressed(const std::string& path) {
    std::map<size_t, size_t> id_to_count;
    for (size_t id : read_to_compressed_) {
        id_to_count[id] ++;
    }

    std::vector<seqan::CharString> compressed_ids;
    std::vector<seqan::Dna5String> compressed_reads;
    size_t current = 0;
    for (size_t i = 0; i < original_ids_.size(); i ++) {
        read_to_compressed_.push_back(current);
        auto id = seqan_string_to_string(original_ids_[i]);
        if (boost::algorithm::ends_with(id, "_copy_1")) {
            current ++;
            compressed_ids.emplace_back((boost::format("cluster___%d___size___%d") % (current - 1) % id_to_count[current - 1]).str());
            compressed_reads.push_back(amplified_reads_[i].read);
        }
    }
    seqan::SeqFileOut compressed_file(path.c_str());
    seqan::writeRecords(compressed_file, compressed_ids, compressed_reads);
}

void PcrSimulator::WriteRcm(const std::string& path, const std::vector<PcrSimulator::Record>& reads, const std::vector<size_t>& map) {
    std::ofstream ofs(path);
    for (size_t i = 0; i < reads.size(); i ++) {
        ofs << reads[i].id;
        if (map[i] != MAP_NOWHERE) {
            ofs << "\t" << map[i];
        }
        ofs << "\n";
    }
}
