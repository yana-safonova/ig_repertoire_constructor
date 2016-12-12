#pragma once

#include "../ig_structs/ig_structure_structs.hpp"

// ----------------------------------------------------------------------------
//  Basic class that adds n-nucleotides in VDJ-recombination
// ----------------------------------------------------------------------------

template<class VDJ_Recombination_Ptr, class NInsertionStrategy, class NInsertionSettings>
class NNucleotidesCreator {
    NInsertionStrategy n_insertion_strategy_;

public:
    NNucleotidesCreator(NInsertionStrategy n_insertion_strategy) :
        n_insertion_strategy_(n_insertion_strategy) { }

    VDJ_Recombination_Ptr CreateNNucleotides(VDJ_Recombination_Ptr vdj_recombination) {
        NInsertionSettings settings = n_insertion_strategy_.CreateNInsertionSettings(vdj_recombination);
        vdj_recombination->AddNInsertionSettings(settings);
        return vdj_recombination;
    }
};

// ----------------------------------------------------------------------------
//      Simple N nucleotides creator
//      - generates n random nucleotides
//      - length of N insertion is uniformly distribution from 0 to MAX_LEN
// ----------------------------------------------------------------------------

// HC
struct HC_SimpleNInsertionConstants {
    const static size_t vd_insertion_max_len = 10;
    const static size_t dj_insertion_max_len = 10;
};

class HC_SimpleNInsertionStrategy {
public:
    HC_NInsertionSettings CreateNInsertionSettings(HC_VDJ_Recombination_Ptr vdj_recombination) {
        size_t vd_insertion_len = RandomInt(HC_SimpleNInsertionConstants::vd_insertion_max_len);
        size_t dj_insertion_len = RandomInt(HC_SimpleNInsertionConstants::dj_insertion_max_len);
        return HC_NInsertionSettings(GetRandomSequence(vd_insertion_len), GetRandomSequence(dj_insertion_len));
    }
};

// LC
struct LC_SimpleNInsertionConstants {
    const static size_t vj_insertion_max_len = 10;
};

class LC_SimpleNInsertionStrategy {
public:
    LC_NInsertionSettings CreateNInsertionSettings(LC_VDJ_Recombination_Ptr vdj_recombination) {
        size_t vj_insertion_len = RandomInt(LC_SimpleNInsertionConstants::vj_insertion_max_len);
        return LC_NInsertionSettings(GetRandomSequence(vj_insertion_len));
    }
};

// ----------------------------------------------------------------------------

typedef NNucleotidesCreator<HC_VDJ_Recombination_Ptr, HC_SimpleNInsertionStrategy, HC_NInsertionSettings>
    HC_SimpleNNucleotidesCreator;

typedef NNucleotidesCreator<LC_VDJ_Recombination_Ptr, LC_SimpleNInsertionStrategy, LC_NInsertionSettings>
        LC_SimpleNNucleotidesCreator;