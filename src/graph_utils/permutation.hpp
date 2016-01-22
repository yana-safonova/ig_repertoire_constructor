#pragma once

#include "include_me.hpp"

class Permutation {
    size_t num_vertices_;
    vector<size_t> direct_;
    vector<size_t> reverse_;

public:
    Permutation(size_t num_vertices) :
        num_vertices_(num_vertices) {
        for(size_t i = 0; i < num_vertices_; i++) {
            direct_.push_back(i);
            reverse_.push_back(i);
        }
    }

    void ReadFromFile(string filename);

    size_t Size() const { return num_vertices_; }

    const vector<size_t>& Direct() const { return direct_; }

    const vector<size_t>& Reverse() const { return reverse_; }

private:
    DECL_LOGGER("Permutation");
};

ostream& operator<<(ostream &out, const Permutation &permutation);

typedef shared_ptr<Permutation> PermutationPtr;
