//
// Created by Andrew Bzikadze on 4/12/17.
//

#pragma once

#include <memory>
#include "ig_simulator_utils.hpp"
#include "clonal_trees/tree/node.hpp"
#include <random>

namespace ig_simulator {

class AbstractShmCreator {
public:
    AbstractShmCreator() = default;
    AbstractShmCreator(const AbstractShmCreator&) = delete;
    AbstractShmCreator(AbstractShmCreator&&) = delete;
    AbstractShmCreator& operator=(const AbstractShmCreator&) = delete;
    AbstractShmCreator& operator=(AbstractShmCreator&&) = delete;

    virtual ~AbstractShmCreator() { }

    virtual Node::SHM_Vector GenerateSHM_Vector(size_t length) const = 0;
};

using AbstractShmCreatorCPtr = std::unique_ptr<AbstractShmCreator>;


class PoissonShmCreator final : public AbstractShmCreator {
private:
    mutable std::poisson_distribution<size_t> distribution;

public:
    PoissonShmCreator(double lambda):
        AbstractShmCreator(),
        distribution(check_numeric_positive(lambda))
    { }


    Node::SHM_Vector GenerateSHM_Vector(size_t length) const override;
};

} // End namespace ig_simulator