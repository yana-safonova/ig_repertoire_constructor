#include "annotated_tree_storage.hpp"

namespace antevolo {
    const AnnotatedEvolutionaryTree& AnnotatedTreeStorage::operator[](size_t index) const {
        VERIFY_MSG(index < size(), "Index " << index << " exceeds storage size: " << size());
        return annotated_trees_[index];
    }
}