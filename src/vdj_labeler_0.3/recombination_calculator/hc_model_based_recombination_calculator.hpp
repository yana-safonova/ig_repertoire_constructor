#pragma once

#include "abstract_recombination_calculator.hpp"
#include "model/recombination_model.hpp"

namespace vdj_labeler {

class HCModelBasedRecombinationCalculator:
        public AbstractRecombinationCalculator<recombination_utils::HCRecombination>
{
    HCProbabilityRecombinationModel model_;

public:
    HCModelBasedRecombinationCalculator(const HCProbabilityRecombinationModel &model) :
        model_(model) { }

    double ComputeAssemblyProbability(const recombination_utils::HCRecombination &recombination) const;
};

} // End namespace vdj_labeler
