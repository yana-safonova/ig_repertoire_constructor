#pragma once

#include "../../ig_structs/ig_structure_structs.hpp"
#include "../../repertoire.hpp"

// RGYM - WRCY based SHM
template<class IgVariableRegionPtr, class SHMSettings>
class RgywWrcySHMStrategy {
    size_t min_mutation_number_;
    size_t max_mutation_number_;
    double substitution_prop_;

    bool IsPosR(const string &seq, size_t pos) {
        return seq[pos] == 'A' or seq[pos] == 'G';
    }

    bool IsPosY(const string &seq, size_t pos) {
        return seq[pos] == 'C' or seq[pos] == 'T';
    }

    bool IsPosW(const string &seq, size_t pos) {
        return seq[pos] == 'A' or seq[pos] == 'T';
    }

    bool IsRGYW(const string &seq, size_t pos) {
        // A - G, G, C - T, A - T
        if(pos == 0 or pos >= seq.size() - 2)
            return false;
        return IsPosR(seq, pos - 1) and seq[pos] == 'G' and IsPosY(seq, pos + 1) and IsPosW(seq, pos + 2);
    }

    bool IsWRCY(const string &seq, size_t pos) {
        // A - T, A - G, C, C - T
        if(pos <= 1 or pos == seq.size() -1)
            return false;
        return IsPosW(seq, pos - 2) and IsPosR(seq, pos - 1) and seq[pos] == 'C' and IsPosY(seq, pos + 1);
    }

    bool IsMutationSubstitution() {
        return double(rand()) / RAND_MAX <= substitution_prop_;
    }

    SHM CreateRandomSubstitution(const string &seq, size_t pos) {
        SHM shm(pos, SubstitutionSHM);
        shm.SetSubstitution(GetAnotherRandomNucleotide(seq, pos));
        return shm;
    }

    SHM CreateMutation(const string &seq, size_t pos) {
        if(IsMutationSubstitution())
            return CreateRandomSubstitution(seq, pos);
        SHM shm(pos, DeletionSHM);
        shm.SetNumDeleted(1);
        return shm;
    }

    vector<size_t> HotSpotVector(const string &seq) {
        vector<size_t> hot_spot_pos;
        for(size_t i = 0; i < seq.size(); i++)
            if(IsRGYW(seq, i) or IsWRCY(seq, i)) {
                if(hot_spot_pos.size() == 0)
                    hot_spot_pos.push_back(i);
                else if(i - hot_spot_pos[hot_spot_pos.size() - 1] >= 3)
                    hot_spot_pos.push_back(i);
            }
        return hot_spot_pos;
    }

    size_t GetRandomPosition(const vector<size_t>& hotspot_pos, const set<size_t> &mutated_pos) {
        size_t random_pos = rand() % hotspot_pos.size();
        while(mutated_pos.find(random_pos) != mutated_pos.end())
            random_pos = rand() % hotspot_pos.size();
        return hotspot_pos[random_pos];
    }

    set<size_t> RandomlyAddMutations(const vector<size_t>& pos) {
        set<size_t> mutated_pos_set;
        size_t num_mutations = min<size_t>(pos.size(),
                                           max<size_t>(min_mutation_number_, rand() % max_mutation_number_));
        for(size_t i = 0; i < num_mutations; i++) {
            size_t mutated_pos = GetRandomPosition(pos, mutated_pos_set);
            mutated_pos_set.insert(mutated_pos);
        }
        return mutated_pos_set;
    }

public:
    RgywWrcySHMStrategy(size_t min_mutation_number,
                        size_t max_mutation_number,
                        double substitution_prop = .8) :
            min_mutation_number_(min_mutation_number),
            max_mutation_number_(max_mutation_number),
            substitution_prop_(substitution_prop) { }

    SHMSettings CreateSHM(IgVariableRegionPtr ig_variable_region_ptr) {
        SHMSettings settings;
        string ab_string = ig_variable_region_ptr->Sequence();
        vector<size_t> hotspot_pos = HotSpotVector(ab_string);
        set<size_t> mutation_pos = RandomlyAddMutations(hotspot_pos);
        for(auto it = mutation_pos.begin(); it != mutation_pos.end(); it++)
            settings.Add(CreateMutation(ab_string, *it));
        return settings;
    }
};

typedef RgywWrcySHMStrategy<HC_VariableRegionPtr, SHMSettings> HC_RgywWrcySHMStrategy;
typedef RgywWrcySHMStrategy<LC_VariableRegionPtr, SHMSettings> LC_RgywWrcySHMStrategy;