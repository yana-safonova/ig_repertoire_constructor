#include <verify.hpp>

#include "evolutionary_tree_storage.hpp"

namespace antevolo {
    void EvolutionaryTreeStorage::Add(EvolutionaryTree tree) {
        trees_.push_back(tree); // probably, we can use move here
    }

    void EvolutionaryTreeStorage::AppendArchive(const EvolutionaryTreeStorage& tree_storage) {
        for(auto it = tree_storage.cbegin(); it != tree_storage.cend(); it++) {
            trees_.push_back(*it);
        }
    }

    const EvolutionaryTree& EvolutionaryTreeStorage::GetTreeByIndex(size_t index) const {
        VERIFY_MSG(index < size(), "Index " << index << " exceeds size of evolutionary tree storage: " << size());
        return trees_[index];
    }

    const EvolutionaryTree& EvolutionaryTreeStorage::operator[](size_t index) const {
        return GetTreeByIndex(index);
    }
}