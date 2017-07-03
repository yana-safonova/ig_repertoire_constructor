//
// Created by Andrew Bzikadze on 4/9/17.
//

#pragma once

#include "node.hpp"
#include "base_repertoire/metaroot/metaroot.hpp"

namespace ig_simulator {

class Tree {
    const AbstractMetaroot* metaroot;
    const std::vector<Node> nodes;
    const std::vector<std::string> sequences;

public:
    Tree(const AbstractMetaroot* const metaroot,
         std::vector<Node>&& nodes = {},
         std::vector<std::string>&& sequences = {}) noexcept:
        metaroot(metaroot),
        nodes(std::move(nodes)),
        sequences(std::move(sequences))
    { }

    Tree(const Tree&) = default;
    Tree(Tree&&) = default;

    Tree& operator=(const Tree&) = default;
    Tree& operator=(Tree&&) = default;

    size_t Size() const { return nodes.size(); }
    const AbstractMetaroot* Metaroot() const { return metaroot; }

    const std::vector<std::string>& Sequences() const { return sequences; }

    bool IsNodeIncluded(size_t node_ind) const {
        return nodes[node_ind].IsIncluded();
    }

    friend std::ostream& operator<<(std::ostream& out, const Tree& tree);
};

std::ostream& operator<<(std::ostream& out, const Tree& tree);

} // End namespace ig_simulator
