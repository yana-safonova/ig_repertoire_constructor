//
// Created by Andrew Bzikadze on 4/9/17.
//

#pragma once

#include <cstddef>
#include "cartesian_tree.hpp"
#include "simulation_routines.hpp"
#include "ig_simulator_utils.hpp"

namespace ig_simulator {

class AbstractPoolManager {
protected:
    Treap<> pool;
    std::bernoulli_distribution ret_to_pool_distr;
    size_t max_index;

public:
    AbstractPoolManager(double ret_prob):
            pool(),
            ret_to_pool_distr(check_numeric_positive(ret_prob)),
            max_index(1)
    {
        pool.Insert(0, 1);
        // pool.Insert(1, 1);
    }

    AbstractPoolManager(const AbstractPoolManager&) = delete;
    AbstractPoolManager(AbstractPoolManager&&) = delete;
    AbstractPoolManager& operator=(const AbstractPoolManager&) = delete;
    AbstractPoolManager& operator=(AbstractPoolManager&&) = delete;

    virtual std::pair<size_t, bool> GetIndex() = 0;
};

using AbstractPoolManagerCPtr = std::unique_ptr<AbstractPoolManager>;


class UniformPoolManager final : public AbstractPoolManager {
public:
    UniformPoolManager(double ret_prob):
        AbstractPoolManager(ret_prob)
    { }

    std::pair<size_t, bool> GetIndex() override;
};

class WideTreePoolManager final : public AbstractPoolManager {
public:
    WideTreePoolManager(double ret_prob):
        AbstractPoolManager(ret_prob)
    { }

    std::pair<size_t, bool> GetIndex() override;
};

class DeepTreePoolManager final : public AbstractPoolManager {
public:
    DeepTreePoolManager(double ret_prob):
        AbstractPoolManager(ret_prob)
    { }

    std::pair<size_t, bool> GetIndex() override;
};

} // End namespace ig_simulator
