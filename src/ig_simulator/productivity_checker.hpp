//
// Created by Andrew Bzikadze on 12/21/16.
//

#pragma once

#include "ig_structs/ig_structure_structs.hpp"
#include <cstring>

#include <algorithm>
#include <array>
#include <set>

namespace ig_simulator {

class HC_ProductivityChecker {
private:
    // constexpr static std::array<char[], 3> stop_codons_ {"TAG", "TAA", "TGA"};
    const static size_t hash_base = 10;
    const static set<size_t> stop_codons_hashes_;
    const static std::string nucl_bases;

    size_t StringToHash(const string &s) {
        size_t hash = 0;
        size_t p = 1;
        for (const char nucl : s) {
            hash += nucl_bases.find(nucl) * p;
            p *= hash_base;
        }
        return hash;
    }

    bool IsVJJunctionOutFrame(const HC_VDJ_Recombination_Ptr &vdj) {
        size_t last_V_coding_position = vdj->VgeneLen();
        size_t first_J_coding_position = last_V_coding_position +
                vdj->PInsertionSettings().VEnd().size() + vdj->NInsertionSettings().VD_Insertion().size() +
                vdj->PInsertionSettings().DStart().size() +
                vdj->DgeneLen() +
                vdj->PInsertionSettings().DEnd().size() + vdj->NInsertionSettings().DJ_Insertion().size() +
                vdj->PInsertionSettings().JStart().size();
        std::cout << last_V_coding_position << " " << first_J_coding_position << std::endl;
        return (last_V_coding_position % 3) != (first_J_coding_position % 3);
    }

    bool HasStopCodon(const HC_VDJ_Recombination_Ptr &vdj) {
        const auto &s = vdj->Sequence();
        bool has_stop_codon = false;
        for (size_t i = 0; (not has_stop_codon) and (i + 2 < s.size()); i += 3) {
            size_t substr_hash = StringToHash(s.substr(i, 3));
            has_stop_codon = (stop_codons_hashes_.find(substr_hash) != stop_codons_hashes_.end());
        }
        return has_stop_codon;
    }

public:
    bool IsNonProductive(const HC_VDJ_Recombination_Ptr &vdj) {
        return IsVJJunctionOutFrame(vdj) or HasStopCodon(vdj);
    }
};

const set<size_t> HC_ProductivityChecker::stop_codons_hashes_{203, 3, 23}; // TAG, TAA, TGA
const std::string HC_ProductivityChecker::nucl_bases = "ACGT";

} // End namespace ig_simulator