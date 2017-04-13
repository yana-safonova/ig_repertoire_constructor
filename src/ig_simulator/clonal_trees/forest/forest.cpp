//
// Created by Andrew Bzikadze on 4/9/17.
//

#include "forest.hpp"

namespace ig_simulator {

std::ostream& operator<<(std::ostream& out, const Forest& forest) {
    for(size_t i = 0; i < forest.trees.size(); ++i) {
        const auto& tree = forest.trees[i];
        out << "===============================================\n";
        out << "Tree    # " << i + 1 << " / " << forest.trees.size() << '\n';
        out << "===============================================\n";
        out << tree;
        out << '\n';
    }
    return out;
}

} // End namespace ig_simulator
