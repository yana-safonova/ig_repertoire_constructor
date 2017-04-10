//
// Created by Andrew Bzikadze on 4/9/17.
//

#pragma once

#include <vector>
#include <utility>
#include <seqan/basic.h>

namespace ig_simulator {

class Node {
public:
    using SHM_Vector = std::vector<std::pair<size_t, seqan::Dna>>;

private:
    const Node* parent;

    // We store only SHMs "on the edge" from the parent
    SHM_Vector shms;


public:
    Node(const Node* parent = nullptr, SHM_Vector&& shms = {}):
        parent(parent),
        shms(shms)
    { }

    Node(const Node&) = default;

    // I redefine move-constructor as I want node.parent to be nullptr afterwards.
    Node(Node&& node) noexcept:
        parent(std::move(node.parent)),
        shms(std::move(node.shms))
    {
        node.parent = nullptr;
    }

    Node& operator=(const Node& node) = default;

    // I redefine move-assignment operator as I want node.parent to be nullptr afterwards.
    Node& operator=(Node&& node) {
        swap(node);
        node.parent = nullptr;
        return *this;
    }

    void swap(Node& node) {
        std::swap(parent, node.parent);
        std::swap(shms, node.shms);
    }
};

} // End namespace ig_simulator
