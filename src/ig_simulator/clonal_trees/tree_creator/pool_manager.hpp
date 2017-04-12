//
// Created by Andrew Bzikadze on 4/9/17.
//

#pragma once

#include <cstddef>
#include "cartesian_tree.hpp"
#include "random_index.hpp"

namespace ig_simulator {

class AbstractPoolManager {
protected:
    Treap<> pool;

public:
    AbstractPoolManager(): pool() {
        pool.Insert(0, 1);
    }

    AbstractPoolManager(const AbstractPoolManager&) = delete;
    AbstractPoolManager(AbstractPoolManager&&) = delete;
    AbstractPoolManager& operator=(const AbstractPoolManager&) = delete;
    AbstractPoolManager& operator=(AbstractPoolManager&&) = delete;

    virtual size_t GetIndex() = 0;
};

class UniformPoolManager final : public AbstractPoolManager {
public:
    size_t GetIndex() override {
        size_t sum = pool.Sum();
        size_t raw_index = random_index(1, sum);
        size_t index = pool.Find(raw_index);
        pool.Insert(pool.Size(), 1);
        // pool.Erase(index, 1);
        return index;
    }
};

} // End namespace ig_simulator
