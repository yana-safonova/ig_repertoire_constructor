#pragma once

#include "../model/recombination_model.hpp"
#include "../recombination/recombination.hpp"

template<class Recombination>
class AbstractRecombinationCalculator {
public:
    virtual double ComputeAssemblyProbability(Recombination recombination) = 0;
    virtual ~AbstractRecombinationCalculator() { }
};