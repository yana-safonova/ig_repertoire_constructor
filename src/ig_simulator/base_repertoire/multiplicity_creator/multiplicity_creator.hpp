//
// Created by Andrew Bzikadze on 3/27/17.
//

#pragma once

#include <cstddef>
#include <random>

#include "verify.hpp"

#include "random_generator.hpp"
#include "ig_simulator_config.hpp"

namespace ig_simulator {

class AbstractMultiplicityCreator {
public:
    virtual size_t RandomMultiplicity() = 0;
};

using AbstractMultiplicityCreatorPtr = std::unique_ptr<AbstractMultiplicityCreator>;

class GeometricMultiplicityCreator : public AbstractMultiplicityCreator {
private:
    double lambda;
    std::geometric_distribution<size_t> distribution;

    static double check_lambda(double lambda) {
        VERIFY(lambda > 0);
        return lambda;
    }

public:
    GeometricMultiplicityCreator(double lambda):
        lambda(lambda),
        distribution(check_lambda(lambda))
    { }

    GeometricMultiplicityCreator(const MultiplicityCreatorParams::GeometricParams &config):
        GeometricMultiplicityCreator(config.lambda)
    { }

    virtual size_t RandomMultiplicity() override {
        return distribution(MTSingleton::GetInstance()) + 1;
    }

    double Mean() const { return 1. / lambda + 1; }
};

AbstractMultiplicityCreatorPtr get_multiplicity_creator(const MultiplicityCreatorParams &config);

} // End namespace ig_simulator
