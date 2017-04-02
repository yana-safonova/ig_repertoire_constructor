//
// Created by Andrew Bzikadze on 4/2/17.
//

#pragma once

#include "metaroot/metaroot.hpp"

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
};

} // End namespace ig_simulator
