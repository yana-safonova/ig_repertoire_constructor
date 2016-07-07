#pragma once

#include "model/recombination_model.hpp"
#include "recombination_utils/recombination.hpp"

namespace vdj_labeler {

template<class Recombination>
class AbstractRecombinationCalculator {
public:
    virtual double ComputeAssemblyProbability(const Recombination &recombination) const = 0;
    virtual ~AbstractRecombinationCalculator() { }
};

} // End namespace vdj_labeler
