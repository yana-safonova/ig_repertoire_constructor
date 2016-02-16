#include "logger/logger.hpp"
#include "raw_pairing_data_storage.hpp"
#include "pairing_fastq_utils.hpp"

void IsotypeUmiSequences::CheckConsistency(const UmiIsotypeSequence &umi_sequence) const {
    assert(isotype_ == umi_sequence.isotype);
}

void IsotypeUmiSequences::Update(UmiIsotypeSequence umi_sequence) {
    CheckConsistency(umi_sequence);
    sequences_.push_back(umi_sequence);
}

void RawPairingDataStorage::UpdateRecord(std::string db, std::string header, std::string sequence) {
    if(db == "" or header == "" or sequence == "")
        return;
    if(db_index_map_.find(db) == db_index_map_.end()) {
        RawPairingDataPtr pairing_record(new RawPairingData(db));
        raw_pairing_records_.push_back(pairing_record);
        db_index_map_[db] = raw_pairing_records_.size() - 1;
    }
    size_t index = db_index_map_[db];
    raw_pairing_records_[index]->Update(header, sequence);
}

void RawPairingDataStorage::Update(std::string fastq_fname) {
    INFO("Adding pairing data from " << fastq_fname);
    size_t num_lines_before = raw_pairing_records_.size();
    std::ifstream fastq_fhandler(fastq_fname);
    assert(fastq_fhandler.good());
    std::string db;
    std::string header;
    std::string sequence;
    bool plus_was_read = false;
    while(!fastq_fhandler.eof()) {
        std::string tmp;
        getline(fastq_fhandler, tmp);
        if(PairingFastqUtils::LineIsHeader(tmp)) {
            UpdateRecord(db, header, sequence);
            db = PairingFastqUtils::ExtractDropletBarcode(tmp);
            header = tmp;
            sequence = "";
            plus_was_read = false;
        }
        else if(PairingFastqUtils::LineIsPlus(tmp))
            plus_was_read = true;
        else if(!plus_was_read)
            sequence += tmp;
    }
    size_t num_lines_after = raw_pairing_records_.size();
    INFO(num_lines_after - num_lines_before << " new records were added");
}

RawPairingDataPtr RawPairingDataStorage::operator[](size_t index) {
    assert(index < raw_pairing_records_.size());
    return raw_pairing_records_[index];
}

RawPairingDataPtr RawPairingDataStorage::GetRecordByDb(std::string db) {
    assert(db_index_map_.find(db) != db_index_map_.end());
    return raw_pairing_records_[db_index_map_[db]];
}