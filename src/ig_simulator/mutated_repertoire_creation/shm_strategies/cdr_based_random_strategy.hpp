#pragma once

#include "../../ig_structs/ig_structure_structs.hpp"
#include "../../repertoire.hpp"

template<class IgVariableRegionPtr, class SHMSettings>
class CDRBasedRandomSHMStrategy {
    size_t min_mutation_number_;
    size_t max_mutation_number_;
    double fr_mutation_prop_;

    size_t GetRandomMutationNumber() {
        return RandomIndex(max_mutation_number_, min_mutation_number_);
    }

    bool IsMutationInCDR() {
        return double(rand()) / RAND_MAX > fr_mutation_prop_;
    }

    SHM CreateMutationInCDR(IgVariableRegionPtr ig_variable_region_ptr) {
        size_t cdr_index = RandomInt(3);
        size_t mutation_index = size_t(-1);
        if(cdr_index == 0) // cdr1
            mutation_index = RandomIndex(ig_variable_region_ptr->GetCDRSettings().CDR1().End(),
                                         ig_variable_region_ptr->GetCDRSettings().CDR1().Start());
        else if(cdr_index == 1) // cdr2
            mutation_index = RandomIndex(ig_variable_region_ptr->GetCDRSettings().CDR2().End(),
                                         ig_variable_region_ptr->GetCDRSettings().CDR2().Start());
        else // cdr3
            mutation_index = RandomIndex(ig_variable_region_ptr->GetCDRSettings().CDR3().End(),
                                         ig_variable_region_ptr->GetCDRSettings().CDR3().Start());
        SHM mutation(mutation_index, SubstitutionSHM);
        mutation.SetSubstitution(GetAnotherRandomNucleotide(ig_variable_region_ptr->Sequence()[mutation_index]));
        return mutation;
    }

    SHM CreateMutationInFR(IgVariableRegionPtr ig_variable_region_ptr) {
        size_t fr_index = RandomInt(4); // index of one of 4 frs
        size_t mutation_index = size_t(-1);
        if(fr_index == 0)
            mutation_index = RandomIndex(ig_variable_region_ptr->GetCDRSettings().CDR1().Start());
        else if(fr_index == 1)
            mutation_index = RandomIndex(ig_variable_region_ptr->GetCDRSettings().CDR2().Start(),
                                         ig_variable_region_ptr->GetCDRSettings().CDR1().End());
        else if(fr_index == 2)
            mutation_index = RandomIndex(ig_variable_region_ptr->GetCDRSettings().CDR3().Start(),
                                         ig_variable_region_ptr->GetCDRSettings().CDR2().End());
        else if(fr_index == 3)
            mutation_index = RandomIndex(ig_variable_region_ptr->Length(),
                                         ig_variable_region_ptr->GetCDRSettings().CDR3().End());
        SHM mutation(mutation_index, SubstitutionSHM);
        mutation.SetSubstitution(GetAnotherRandomNucleotide(ig_variable_region_ptr->Sequence()[mutation_index]));
        return mutation;
    }

public:
    CDRBasedRandomSHMStrategy(size_t min_mutation_number, size_t max_mutation_number, double fr_mutation_prop = .25) :
        min_mutation_number_(min_mutation_number),
        max_mutation_number_(max_mutation_number),
        fr_mutation_prop_(fr_mutation_prop) { }

    SHMSettings CreateSHM(IgVariableRegionPtr ig_variable_region_ptr) {
        SHMSettings shm_settings;
        string ab_string = ig_variable_region_ptr->Sequence();
        size_t num_mutations = GetRandomMutationNumber();
        for (size_t i = 0; i < num_mutations; i++)
            if (IsMutationInCDR()) {
                shm_settings.Add(CreateMutationInCDR(ig_variable_region_ptr));
            }
            else {
                shm_settings.Add(CreateMutationInFR(ig_variable_region_ptr));
            }
        return shm_settings;
    }
};

typedef CDRBasedRandomSHMStrategy<HC_VariableRegionPtr, SHMSettings> HC_CDRBasedRandomSHMStrategy;
typedef CDRBasedRandomSHMStrategy<LC_VariableRegionPtr, SHMSettings> LC_CDRBasedRandomSHMStrategy;