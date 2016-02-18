#pragma once

#include <unordered_set>
#include "../pairing_ig_data/raw_pairing_data_storage.hpp"

struct PairingDataCdr3 {
    std::unordered_set<std::string> hcdr3s;
    std::unordered_set<std::string> kcdr3s;
    std::unordered_set<std::string> lcdr3s;
};

class Cdr3GraphComputer {
    const RawPairingDataStorage &pairing_data_storage_;

    std::unordered_map<std::string, std::vector<size_t>> hcdr3_map_;
    std::unordered_map<std::string, std::vector<size_t>> kcdr3_map_;
    std::unordered_map<std::string, std::vector<size_t>> lcdr3_map_;

    std::vector<PairingDataCdr3> pairing_data_cdr3s_;

    std::vector<std::unordered_set<std::string>> related_hcdr3_;
    std::vector<std::unordered_set<std::string>> related_kcdr3_;
    std::vector<std::unordered_set<std::string>> related_lcdr3_;

    std::vector<size_t> record_processed_;
    std::unordered_map<std::string, bool> hcdr3_processed_;
    std::unordered_map<std::string, bool> kcdr3_processed_;
    std::unordered_map<std::string, bool> lcdr3_processed_;

    void Initialize();

    void ComputeHCdr3s();

    void ComputeKCdr3s();

    void ComputeLCdr3s();

    void InitializeMaps();

    std::string GetStartHCdr3();

    std::string GetStartKCdr3();

    std::string GetStartLCdr3();

    void ComputeCdr3GraphFromStartCdr3(std::string hcdr3);

    void ComputeCdr3Graphs();

public:
    Cdr3GraphComputer(const RawPairingDataStorage &pairing_data_storage) :
            pairing_data_storage_(pairing_data_storage) {
        Initialize();
    }

    void Compute();
};