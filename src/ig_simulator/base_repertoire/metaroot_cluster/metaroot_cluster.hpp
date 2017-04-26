//
// Created by Andrew Bzikadze on 4/2/17.
//

#pragma once

#include "base_repertoire/metaroot/metaroot.hpp"

namespace ig_simulator {

class MetarootCluster {
private:
    AbstractMetarootCPtr metaroot_p;
    size_t multiplicity;

public:
    MetarootCluster(AbstractMetarootCPtr&& metaroot_p,
                    const size_t multiplicity):
        metaroot_p(std::move(metaroot_p)),
        multiplicity(multiplicity)
    { }

    MetarootCluster(const MetarootCluster&) = delete;
    MetarootCluster(MetarootCluster&&) = default;
    MetarootCluster& operator=(const MetarootCluster&) = delete;
    MetarootCluster& operator=(MetarootCluster&&) = delete;

    const AbstractMetarootCPtr& MetarootPtr() const { return metaroot_p; }
    size_t Multiplicity() const { return multiplicity; }
};

} // End namespace ig_simulator
