//
// Created by Andrew Bzikadze on 3/27/17.
//

#pragma once

#include <cstddef>
#include "random_generator.hpp"
#include <random>
#include "verify.hpp"

namespace ig_simulator {

class AbstractMultiplicityCreator {
public:
    virtual size_t RandomMultiplicity() const = 0;
};

class GeometricMultiplicityCreator : public AbstractMultiplicityCreator {
private:
    double lambda;
    std::exponential_distribution<size_t> distribution;

public:
    GeometricMultiplicityCreator(double lambda):
        lambda(lambda),
        distribution(lambda)
    {
        VERIFY(lambda > 0);
    }

    virtual size_t RandomMultiplicity() const override {
        return distribution(MTSingleton::GetInstance()) + 1;
    }

    double Mean() const { return 1. / lambda + 1; }

};

} // End namespace ig_simulator
