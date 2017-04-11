//
// Created by Andrew Bzikadze on 4/9/17.
//

#pragma once

#include <cstddef>
#include <vector>

namespace ig_simulator {

class AbstractPoolManager {
    std::vector<size_t> pool;

public:
    AbstractPoolManager() = default;
    AbstractPoolManager(const AbstractPoolManager&) = delete;
    AbstractPoolManager(AbstractPoolManager&&) = delete;
    AbstractPoolManager& operator=(const AbstractPoolManager&) = delete;
    AbstractPoolManager& operator=(AbstractPoolManager&&) = delete;
};

} // End namespace ig_simulator
