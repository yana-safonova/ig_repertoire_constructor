#pragma once

#include "abstract_recombination_calculator.hpp"

class HCModelBasedRecombinationCalculator : public AbstractRecombinationCalculator<HCRecombination> {
    HCProbabilityRecombinationModel model_;

public:
    HCModelBasedRecombinationCalculator(HCProbabilityRecombinationModel model) :
            model_(model) { }

    double ComputeAssemblyProbability(HCRecombination recombination);
};