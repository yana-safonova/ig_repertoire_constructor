#pragma once

#include <random>
#include <seqan/file.h>
#include "options.hpp"

class PcrSimulator {
private:
    static std::default_random_engine random_engine_;

    struct Record {
        seqan::CharString id;
        const seqan::Dna5String read;
        const seqan::Dna5String barcode;
        const size_t error_count;

        Record(const seqan::CharString& id_, const seqan::Dna5String& read_, const seqan::Dna5String& barcode_, size_t error_count_) :
                id(id_), read(read_), barcode(barcode_), error_count(error_count_) {}
    };

    const Options::SimulationOptions options_;
    std::vector<seqan::CharString> original_ids_;
    std::vector<seqan::Dna5String> original_reads_;
    std::vector<Record> amplified_reads_;
    std::vector<size_t> perm_;
    std::vector<size_t> read_to_compressed_;
    std::vector<size_t> read_to_original_;
    size_t barcode_error_count_;

public:
    static const size_t CHIMERIC_READ_ERROR_COUNT;

    PcrSimulator(const Options::SimulationOptions& options) : options_(options) {
        random_engine_.seed(8356);
    };
    void ReadRepertoire(const std::string& repertoire_file_path);
    void Amplify(size_t output_estimation_limit);
    void WriteResults(const std::string& output_dir_path);

private:
    static const size_t MAP_NOWHERE;

    void CheckLimit(size_t output_estimation_limit);
    seqan::Dna5String GenerateBarcode();
    void ReportAverageErrorRate();
    void SimulatePcr();
    void AddRecord(const std::string& id, const seqan::Dna5String& read, const seqan::Dna5String& barcode, size_t error_count);
    void AddChimericRecord(const std::string& id, const seqan::Dna5String& read, const seqan::Dna5String& barcode);
    void AmplifySequences(double pcr_error_prob);
    void UpdateReadIds();
    void WriteRepertoire(const std::string& path);
    void WriteCompressed(const std::string& path);
    void WriteRcm(const std::string& path, const std::vector<Record>& reads, const std::vector<size_t>& map);
};
