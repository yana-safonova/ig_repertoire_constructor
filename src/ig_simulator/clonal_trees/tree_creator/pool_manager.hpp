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

public:
    AbstractPoolManager(double ret_prob):
            pool(),
            ret_to_pool_distr(check_numeric_positive(ret_prob))
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

class UniformPoolManager final : public AbstractPoolManager {
public:
    UniformPoolManager(double ret_prob):
        AbstractPoolManager(ret_prob)
    { }

    std::pair<size_t, bool> GetIndex() override {
        double raw_index = uniform_double(0, pool.Sum());
        size_t index, freq;
        std::tie(index, freq) = pool.LowerBound(raw_index);
        VERIFY(freq == 1);
        pool.Insert(pool.Size(), 1);
        bool ret_to_pool = ret_to_pool_distr(MTSingleton::GetInstance());
        if (not ret_to_pool) {
            pool.Erase(index, freq);
        }
        return { index, ret_to_pool };
    }
};

class WideTreePoolManager final : public AbstractPoolManager {
public:
    WideTreePoolManager(double ret_prob):
        AbstractPoolManager(ret_prob)
    { }

    std::pair<size_t, bool> GetIndex() override {
        double raw_index = uniform_double(0, pool.Sum());
        size_t index, freq;
        std::tie(index, freq) = pool.LowerBound(raw_index);
        pool.Insert(pool.Size(), 1);
        bool ret_to_pool = ret_to_pool_distr(MTSingleton::GetInstance());
        pool.Erase(index, freq);
        if (ret_to_pool) {
            pool.Insert(index, freq + 1);
        }
        return { index, ret_to_pool };
    }
};

class DeepTreePoolManager final : public AbstractPoolManager {
public:
    DeepTreePoolManager(double ret_prob):
        AbstractPoolManager(ret_prob)
    { }

    std::pair<size_t, bool> GetIndex() override {
        double raw_index = uniform_double(0, pool.Sum());
        size_t index, freq;
        std::tie(index, freq) = pool.LowerBound(raw_index);
        pool.Insert(pool.Size(), freq);
        bool ret_to_pool = ret_to_pool_distr(MTSingleton::GetInstance());
        pool.Erase(index, freq);
        if (ret_to_pool) {
            pool.Insert(index, freq + 1);
        }
        return { index, ret_to_pool };
    }
};

} // End namespace ig_simulator
