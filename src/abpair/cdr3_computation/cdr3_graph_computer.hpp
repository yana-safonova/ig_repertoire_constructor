#pragma once

#include <unordered_set>
#include <abpair_config.hpp>
#include "../pairing_ig_data/raw_pairing_data_storage.hpp"

struct Annotation {
    std::string cdr3;
    std::string v;
    std::string j;

    Annotation() :
            cdr3(),
            v(),
            j() {}

    Annotation(const Annotation &obj) :
            cdr3(obj.cdr3),
            v(obj.v),
            j(obj.j) {}

    bool operator ==(const Annotation &obj) const { return cdr3 == obj.cdr3 and v == obj.v and j == obj.j; }
};

struct AnnotationHasher {
    size_t operator()(const Annotation &obj) const;
};

struct PairingDataCdr3 {
    std::unordered_set<Annotation, AnnotationHasher> hcdr3s;
    std::unordered_set<Annotation, AnnotationHasher> kcdr3s;
    std::unordered_set<Annotation, AnnotationHasher> lcdr3s;
};

class Cdr3GraphComputer {
    const abpair_config::io_config &io_config_;
    const RawPairingDataStorage &pairing_data_storage_;

    std::unordered_map<Annotation, std::vector<size_t>, AnnotationHasher> hcdr3_map_;
    std::unordered_map<Annotation, std::vector<size_t>, AnnotationHasher> kcdr3_map_;
    std::unordered_map<Annotation, std::vector<size_t>, AnnotationHasher> lcdr3_map_;

    // cdr3 structure for complete records
    std::vector<PairingDataCdr3> pairing_data_cdr3s_;
    // vector that stores indices in pairing data storage for each complete record
    std::vector<size_t> complete_indices_;
    // vector of collections of related complete records
    std::vector<std::unordered_set<size_t>> related_records_;
    // vector that stores flags "complete record is processed"
    std::vector<bool> record_processed_;

    std::unordered_map<Annotation, bool, AnnotationHasher> hcdr3_processed_;
    std::unordered_map<Annotation, bool, AnnotationHasher> kcdr3_processed_;
    std::unordered_map<Annotation, bool, AnnotationHasher> lcdr3_processed_;

    void Initialize();

    void ComputeHCdr3s();

    void ComputeKCdr3s();

    void ComputeLCdr3s();

    void InitializeMaps();

    size_t GetStartRecordIndex();

    std::unordered_set<size_t> GetRecordsThatShareCdr3sWith(size_t record_ind);

    void ComputeCdr3GraphFromStartRecord(size_t start_record_ind);

    void ComputeCdr3Graphs();

    // the first number from the pair - # HC isotypes
    // the second number from the pair - # LC isotypes
    std::pair<size_t, size_t> ComputeIsotypesNumberForRelatedGroup(size_t index);

    // the first number from the pair - # CDR3 in KC
    // the second number from the pair - # CDR3 in LC
    std::pair<size_t, size_t> ComputeNumberLightCdr3InRelatedGroup(size_t index);

    bool OutputRelatedGroup(size_t index);

    std::string GetRelatedGroupDirName(size_t index, std::pair<size_t, size_t> num_isotypes);

    void OutputRelatedGroupHcs(size_t index, std::string output_dir);

    void OutputRelatedGroupLcs(size_t index, std::string output_dir);

    void OutputRelatedGroupKcs(size_t index, std::string output_dir);

    void OutputGraph(size_t index, std::string output_dir);

public:
    Cdr3GraphComputer(const abpair_config::io_config &io_config,
            const RawPairingDataStorage &pairing_data_storage) :
            io_config_(io_config),
            pairing_data_storage_(pairing_data_storage) {
        Initialize();
    }

    void Compute();

    void OutputCdr3Graphs();
};