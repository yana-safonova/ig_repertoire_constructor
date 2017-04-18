//
// Created by Andrew Bzikadze on 4/9/17.
//

#pragma once

#include "node.hpp"
#include "base_repertoire/metaroot/metaroot.hpp"

namespace ig_simulator {

class Tree {
    const AbstractMetaroot* metaroot;
    std::vector<Node> nodes;

public:
    Tree(const AbstractMetaroot* const metaroot,
         std::vector<Node>&& nodes = {}) noexcept:
        metaroot(metaroot),
        nodes(std::move(nodes))
    { }

    Tree(const Tree&) = default;
    Tree(Tree&&) = default;

    Tree& operator=(const Tree&) = default;
    Tree& operator=(Tree&&) = default;

    size_t Size() const { return nodes.size(); }
    const AbstractMetaroot* Metaroot() const { return metaroot; }

    std::vector<std::string> Sequences() const {
        std::vector<std::string> sequences;
        sequences.reserve(Size());
        sequences.emplace_back(metaroot->Sequence());
        for(size_t i = 1; i < Size(); ++i) {
            const auto& node = nodes[i];
            sequences.emplace_back(sequences[node.ParentInd()]);

            std::string& seq = sequences.back();
            for(const auto& shm : node.SHMs()) {
                size_t value = seqan::Dna5(seq[shm.first]).value;
                seq[shm.first] = seqan::Dna(value > shm.second ? shm.second : ((shm.second + 1) & 3));
            }
        }

        return sequences;
    }

    bool IsNodeIncluded(size_t node_ind) const {
        return nodes[node_ind].IsIncluded();
    }

    friend std::ostream& operator<<(std::ostream& out, const Tree& tree);
};

std::ostream& operator<<(std::ostream& out, const Tree& tree);

} // End namespace ig_simulator
