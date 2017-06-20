//
// Created by Andrew Bzikadze on 3/27/17.
//

#pragma once

#include <cstddef>
#include <random>

#include "verify.hpp"

#include "ig_simulator_utils.hpp"
#include "random_generator.hpp"
#include "ig_simulator_config.hpp"

namespace ig_simulator {

class AbstractMultiplicityCreator {
public:
    AbstractMultiplicityCreator() = default;
    AbstractMultiplicityCreator(const AbstractMultiplicityCreator&) = delete;
    AbstractMultiplicityCreator(AbstractMultiplicityCreator&&) = delete;
    AbstractMultiplicityCreator& operator=(const AbstractMultiplicityCreator&) = delete;
    AbstractMultiplicityCreator& operator=(AbstractMultiplicityCreator&&) = delete;

    virtual size_t RandomMultiplicity() = 0;
};

using AbstractMultiplicityCreatorPtr = std::unique_ptr<AbstractMultiplicityCreator>;

class GeometricMultiplicityCreator final : public AbstractMultiplicityCreator {
private:
    double lambda;
    std::geometric_distribution<size_t> distribution;

public:
    GeometricMultiplicityCreator(double lambda):
        lambda(lambda),
        distribution(check_numeric_positive(lambda))
    { }

    GeometricMultiplicityCreator(const MultiplicityCreatorParams::GeometricParams &config):
        GeometricMultiplicityCreator(config.lambda)
    { }

    size_t RandomMultiplicity() override {
        return distribution(MTSingleton::GetInstance()) + 1;
    }

    double Mean() const { return 1. / lambda + 1; }
};

AbstractMultiplicityCreatorPtr get_multiplicity_creator(const MultiplicityCreatorParams &config);

} // End namespace ig_simulator
