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
    mutable std::bernoulli_distribution ret_to_pool_distr;
    size_t max_index;

public:
    AbstractPoolManager(double ret_prob):
            pool(),
            ret_to_pool_distr(check_probability(ret_prob)),
            max_index(1)
    {
        pool.Insert(0, 1);
    }

    AbstractPoolManager(const AbstractPoolManager&) = delete;
    AbstractPoolManager(AbstractPoolManager&&) = delete;
    AbstractPoolManager& operator=(const AbstractPoolManager&) = delete;
    AbstractPoolManager& operator=(AbstractPoolManager&&) = delete;

    size_t MaxIndex() const { return max_index; }
    void Erase(size_t index) {
        VERIFY(index < max_index);
        pool.Erase(index);
    }

    size_t Size() const { return pool.Size(); }
    virtual std::pair<size_t, bool> GetIndex(size_t n_insert) = 0;
};

using AbstractPoolManagerCPtr = std::unique_ptr<AbstractPoolManager>;


class UniformPoolManager final : public AbstractPoolManager {
public:
    UniformPoolManager(double ret_prob):
        AbstractPoolManager(ret_prob)
    { }

    std::pair<size_t, bool> GetIndex(size_t n_insert) override;
};

class WideTreePoolManager final : public AbstractPoolManager {
public:
    WideTreePoolManager(double ret_prob):
        AbstractPoolManager(ret_prob)
    { }

    std::pair<size_t, bool> GetIndex(size_t n_insert) override;
};

class DeepTreePoolManager final : public AbstractPoolManager {
public:
    DeepTreePoolManager(double ret_prob):
        AbstractPoolManager(ret_prob)
    { }

    std::pair<size_t, bool> GetIndex(size_t n_insert) override;
};

} // End namespace ig_simulator
