#pragma once

#include "../ig_structs/ig_structure_structs.hpp"

template<class GeneDB_Ptr, class VDJ_Recombination_Ptr, class RecombinationCreator>
class VDJ_Recombinator {
    const GeneDB_Ptr gene_db_ptr_;
    RecombinationCreator recombination_creator_;

public:
    VDJ_Recombinator(const GeneDB_Ptr gene_db_ptr,
            RecombinationCreator recombination_creator) :
            gene_db_ptr_(gene_db_ptr),
            recombination_creator_(recombination_creator) { }

    VDJ_Recombination_Ptr CreateRecombination() {
        return recombination_creator_.CreateRecombination();
    }
};

// ----------------------------------------------------------------------------
//      Simple creator of VDJ recombination
//      - randomly chooses V, D, J genes from database
//      - probabilities from all V, D, and J are equal
// ----------------------------------------------------------------------------

// HC
class HC_SimpleRecombinationCreator {
    const HC_GenesDatabase_Ptr db_ptr_;

public:
    HC_SimpleRecombinationCreator(const HC_GenesDatabase_Ptr db_ptr) :
            db_ptr_(db_ptr) { }

    HC_VDJ_Recombination_Ptr CreateRecombination() {
        return HC_VDJ_Recombination_Ptr(
                new HC_VDJ_Recombination(db_ptr_,
                RandomIndex(db_ptr_->GenesNumber(variable_gene)),
                RandomIndex(db_ptr_->GenesNumber(diversity_gene)),
                RandomIndex(db_ptr_->GenesNumber(join_gene))));
    }
};

// LC
class LC_SimpleRecombinationCreator {
    const LC_GenesDatabase_Ptr db_ptr_;
    size_t seed_;

public:
    LC_SimpleRecombinationCreator(const LC_GenesDatabase_Ptr db_ptr) :
            db_ptr_(db_ptr) { }

    LC_VDJ_Recombination_Ptr CreateRecombination() {
        return LC_VDJ_Recombination_Ptr(
                new LC_VDJ_Recombination(db_ptr_,
                RandomIndex(db_ptr_->GenesNumber(variable_gene)),
                RandomIndex(db_ptr_->GenesNumber(join_gene))));
    }
};

// ----------------------------------------------------------------------------

typedef VDJ_Recombinator<HC_GenesDatabase_Ptr, HC_VDJ_Recombination_Ptr, HC_SimpleRecombinationCreator>
        HC_SimpleRecombinator;
typedef VDJ_Recombinator<LC_GenesDatabase_Ptr, LC_VDJ_Recombination_Ptr, LC_SimpleRecombinationCreator>
        LC_SimpleRecombinator;