#pragma once

#include "../ig_structs/ig_structure_structs.hpp"

template<class VDJ_Recombination_Ptr, class RemovingStrategy, class RemovingSettings>
class ExonucleaseRemover {
    RemovingStrategy remover_;

public:
    ExonucleaseRemover(RemovingStrategy remover) :
        remover_(remover) { }

    VDJ_Recombination_Ptr CreateRemovingSettings(VDJ_Recombination_Ptr vdj_recombination) {
        RemovingSettings removing_settings = remover_.CreateRemovingSettings(vdj_recombination);
        vdj_recombination->AddRemovingSettings(removing_settings);
        return vdj_recombination;
    }
};

// ----------------------------------------------------------------------------
// HC
struct HC_SimpleRemoverConstants {
    const static size_t v_end_max = 10;
    const static size_t d_start_max = 3;
    const static size_t d_end_max = 3;
    const static size_t j_start_max = 4;
};

class HC_SimpleRemovingStrategy {
public:
    HC_SimpleRemovingStrategy() {}

    HC_RemovingSettings CreateRemovingSettings(HC_VDJ_Recombination_Ptr recombination) {
        return HC_RemovingSettings(
                RandomInt(HC_SimpleRemoverConstants::v_end_max),
                RandomInt(HC_SimpleRemoverConstants::d_start_max),
                RandomInt(HC_SimpleRemoverConstants::d_end_max),
                RandomInt(HC_SimpleRemoverConstants::j_start_max)
        );
    }
};

// LC
struct LC_SimpleRemoverConstants {
    const static size_t v_end_max = 10;
    const static size_t j_start_max = 4;
};

class LC_SimpleRemovingStrategy {
public:
    LC_SimpleRemovingStrategy() {}

    LC_RemovingSettings CreateRemovingSettings(LC_VDJ_Recombination_Ptr recombination) {
        return LC_RemovingSettings(
                RandomInt(LC_SimpleRemoverConstants::v_end_max),
                RandomInt(LC_SimpleRemoverConstants::j_start_max)
        );
    }
};
// ----------------------------------------------------------------------------

typedef ExonucleaseRemover<HC_VDJ_Recombination_Ptr, HC_SimpleRemovingStrategy, HC_RemovingSettings>
        HC_SimpleExonucleaseRemover;
typedef ExonucleaseRemover<LC_VDJ_Recombination_Ptr, LC_SimpleRemovingStrategy, LC_RemovingSettings>
        LC_SimpleExonucleaseRemover;