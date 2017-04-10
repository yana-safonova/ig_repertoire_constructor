//
// Created by Andrew Bzikadze on 4/9/17.
//

#pragma once

#include "clonal_trees/tree/tree.hpp"

namespace ig_simulator {

class Forest {
private:
    std::vector<Tree> trees;

public:
    Forest(std::vector<Tree>&& trees = {}) noexcept:
        trees(trees)
    { }

    Forest(const Forest&) = default;
    Forest(Forest&&) = default;

    Forest& operator=(const Forest&) = default;
    Forest& operator=(Forest&&) = default;
};

} // End namespace ig_simulator
