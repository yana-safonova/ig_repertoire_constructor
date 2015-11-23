#pragma once

#include "recombination_generator.hpp"
#include "../recombination/recombination.hpp"

typedef std::vector<HCRecombination>::iterator HCRecombinationIterator;

class BaseHCRecombinationGenerator : public AbstractRecombinationGenerator<HCRecombination, HCRecombinationIterator> {
    std::vector<HCRecombination> recombinations_;

public:
    void ComputeRecombinations();

    HCRecombinationIterator begin() { return recombinations_.begin(); }

    HCRecombinationIterator end() { return recombinations_.end(); }

    // std::size_t size() const { return recombinations_.size(); }
};