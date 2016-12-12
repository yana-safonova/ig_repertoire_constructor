#pragma once

#include "../ig_structs/ig_structure_structs.hpp"

// ----------------------------------------------------------------------------
//  Basic class that adds p-nucleotides in VDJ-recombination
// ----------------------------------------------------------------------------

template<class VDJ_Recombination_Ptr, class PInsertionStrategy, class PInsertionSettings>
class PNucleotidesCreator {
    PInsertionStrategy p_insertion_strategy_;

public:
    PNucleotidesCreator(PInsertionStrategy p_insertion_strategy) :
            p_insertion_strategy_(p_insertion_strategy) { }

    VDJ_Recombination_Ptr CreatePNucleotides(VDJ_Recombination_Ptr vdj_recombination) {
        PInsertionSettings settings = p_insertion_strategy_.CreatePInsertionSettings(vdj_recombination);
        vdj_recombination->AddPInsertionSettings(settings);
        return vdj_recombination;
    }
};

// ----------------------------------------------------------------------------
//      Simple P nucleotides creator
//      - generates length n (half of palyndrome), length is uniformly distributed from 0 to MAX_LEN
//      - generates n random nucleotides
//      - creates palindrome
// ----------------------------------------------------------------------------

// HC
struct HC_SimplePInsertionConstants {
    const static size_t v_end_min = 2;
    const static size_t v_end_max = 4;

    const static size_t d_start_min = 2;
    const static size_t d_start_max = 4;

    const static size_t d_end_min = 2;
    const static size_t d_end_max = 4;

    const static size_t j_start_min = 2;
    const static size_t j_start_max = 4;
};

class HC_SimplePInsertionStrategy {
public:
    HC_PInsertionSettings CreatePInsertionSettings(HC_VDJ_Recombination_Ptr recombination) {
        return HC_PInsertionSettings(
                GetPalindrom(RandomIndex(HC_SimplePInsertionConstants::v_end_max,
                        HC_SimplePInsertionConstants::v_end_min)),
                GetPalindrom(RandomIndex(HC_SimplePInsertionConstants::d_start_max,
                        HC_SimplePInsertionConstants::d_start_min)),
                GetPalindrom(RandomIndex(HC_SimplePInsertionConstants::d_end_max,
                        HC_SimplePInsertionConstants::d_end_min)),
                GetPalindrom(RandomIndex(HC_SimplePInsertionConstants::j_start_max,
                        HC_SimplePInsertionConstants::j_start_min))
        );
    }
};

// LC
struct LC_SimplePInsertionConstants {
    const static size_t v_end_min = 2;
    const static size_t v_end_max = 4;

    const static size_t j_start_min = 2;
    const static size_t j_start_max = 4;
};

class LC_SimplePInsertionStrategy {
public:
    LC_PInsertionSettings CreatePInsertionSettings(LC_VDJ_Recombination_Ptr recombination) {
        return LC_PInsertionSettings(
                GetPalindrom(RandomIndex(LC_SimplePInsertionConstants::v_end_max,
                        LC_SimplePInsertionConstants::v_end_min)),
                GetPalindrom(RandomIndex(LC_SimplePInsertionConstants::j_start_max,
                        LC_SimplePInsertionConstants::j_start_min))
        );
    }
};

// ----------------------------------------------------------------------------

typedef PNucleotidesCreator<HC_VDJ_Recombination_Ptr, HC_SimplePInsertionStrategy, HC_PInsertionSettings>
    HC_SimplePNucleotidesCreator;

typedef PNucleotidesCreator<LC_VDJ_Recombination_Ptr, LC_SimplePInsertionStrategy, LC_PInsertionSettings>
        LC_SimplePNucleotidesCreator;