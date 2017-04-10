//
// Created by Andrew Bzikadze on 4/9/17.
//

#pragma once

#include "node.hpp"

namespace ig_simulator {

class Tree {
    std::vector<Node> nodes;

public:
    Tree(std::vector<Node>&& nodes = {}) noexcept:
        nodes(std::move(nodes))
    { }

    Tree(const Tree&) = default;
    Tree(Tree&&) = default;

    Tree& operator=(const Tree&) = default;
    Tree& operator=(Tree&&) = default;
};

} // End namespace ig_simulator
