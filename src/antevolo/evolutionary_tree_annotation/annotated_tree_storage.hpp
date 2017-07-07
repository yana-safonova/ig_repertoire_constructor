#pragma once

#include "annotated_evolutionary_tree.hpp"

namespace antevolo {
    class AnnotatedTreeStorage {
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;
        std::vector<AnnotatedEvolutionaryTree> annotated_trees_;

    public:
        AnnotatedTreeStorage(const annotation_utils::CDRAnnotatedCloneSet &clone_set) :
                clone_set_(clone_set) { }

        void AddAnnotatedTree(const EvolutionaryTree &tree) {
            annotated_trees_.push_back(AnnotatedEvolutionaryTree(clone_set_, tree));
        }

        const AnnotatedEvolutionaryTree& operator[](size_t index) const;

        typedef std::vector<AnnotatedEvolutionaryTree>::const_iterator AnnotatedTreeIterator;

        AnnotatedTreeIterator cbegin() const { return annotated_trees_.cbegin(); }

        AnnotatedTreeIterator cend() const { return annotated_trees_.cend(); }

        size_t size() const { return annotated_trees_.size(); }
    };
}