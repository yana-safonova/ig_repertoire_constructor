#include "logger/logger.hpp"
#include <verify.hpp>
#include "raw_pairing_data_storage.hpp"
#include "pairing_fastq_utils.hpp"

std::string get_main_gene_name(std::string full_name) {
    std::vector<std::string> splits;
    boost::split(splits, full_name, boost::is_any_of("|"));
    full_name = splits[1];
    size_t pos = full_name.find('*');
    VERIFY_MSG(pos != std::string::npos, "Incorrect gene name: " << full_name);
    return full_name.substr(0, pos);
}

std::pair<std::string, std::string> get_vj_hit_from_line(std::string line) {
    std::vector<std::string> splits;
    boost::split(splits, line, boost::is_any_of("\t"));
    std::string v_gene = get_main_gene_name(splits[4]);
    std::string j_gene = get_main_gene_name(splits[8]);
    //std::cout << v_gene << " " << j_gene << std::endl;
    return std::make_pair(v_gene, j_gene);
}

// temporary
void RawPairingDataStorage::ExtractMap() {
    std::string fname = "SG3M-Fx1_alignment_info.csv";
    std::ifstream vj_fhanlder_(fname);
    std::string line;
    std::getline(vj_fhanlder_, line);
    size_t num = 0;
    size_t threshold = 50000;
    while(!vj_fhanlder_.eof()) {
        std::getline(vj_fhanlder_, line);
        if(line == "")
            break;
        std::vector<std::string> splits;
        boost::split(splits, line, boost::is_any_of("\t"));
        align_record_vj_[splits[0]] = get_vj_hit_from_line(line);
        num++;
        if(num > threshold) {
            INFO(num << " records were extracted");
            threshold += 50000;
        }
    }
    INFO(align_record_vj_.size() << " VJ alignment were extracted from " << fname);
}

void RawPairingDataStorage::UpdateRecord(DropletBarcode db, std::string header, std::string sequence) {
    if(db.IsEmpty() or header == "" or sequence == "")
        return;
    if(db_index_map_.find(db) == db_index_map_.end()) {
        RawPairingDataPtr pairing_record(new RawPairingData(db));
        raw_pairing_records_.push_back(pairing_record);
        db_index_map_[db] = raw_pairing_records_.size() - 1;
    }
    size_t index = db_index_map_[db];

    std::string header2 = header.substr(1, header.size() - 1);
    if(align_record_vj_.find(header2) == align_record_vj_.end())
        return;
    auto vj_hit = align_record_vj_[header2];

    raw_pairing_records_[index]->Update(db, header, sequence, vj_hit.first, vj_hit.second);
}

void RawPairingDataStorage::Update(std::string fastq_fname) {
    INFO("Adding pairing data from " << fastq_fname);
    ExtractMap();
    size_t num_lines_before = raw_pairing_records_.size();
    std::ifstream fastq_fhandler(fastq_fname);
    assert(fastq_fhandler.good());
    std::string db_str;
    std::string header;
    std::string sequence;
    bool plus_was_read = false;
    while(!fastq_fhandler.eof()) {
        std::string tmp;
        getline(fastq_fhandler, tmp);
        if(PairingFastqUtils::LineIsHeader(tmp)) {
            DropletBarcode db = db_factory_.GetDropletBarcodeByFilename(db_str, fastq_fname);
            UpdateRecord(db, header, sequence);
            db_str = PairingFastqUtils::ExtractDropletBarcode(tmp);
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

RawPairingDataPtr RawPairingDataStorage::operator[](size_t index) const {
    assert(index < raw_pairing_records_.size());
    return raw_pairing_records_[index];
}

RawPairingDataPtr RawPairingDataStorage::GetRecordByDb(DropletBarcode db) {
    assert(db_index_map_.find(db) != db_index_map_.end());
    return raw_pairing_records_[db_index_map_[db]];
}