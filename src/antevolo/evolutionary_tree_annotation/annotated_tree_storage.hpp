#pragma once

#include "annotated_evolutionary_tree.hpp"
#include "clone_set_with_fakes.hpp"

namespace antevolo {
    class AnnotatedTreeStorage {
//        const CloneSetWithFakes &clone_set_;
        std::vector<AnnotatedEvolutionaryTree> annotated_trees_;

    public:
        void AddAnnotatedTree(const EvolutionaryTree &tree) {
            annotated_trees_.push_back(AnnotatedEvolutionaryTree(tree));
        }

        const AnnotatedEvolutionaryTree& operator[](size_t index) const;

        typedef std::vector<AnnotatedEvolutionaryTree>::const_iterator AnnotatedTreeIterator;

        AnnotatedTreeIterator cbegin() const { return annotated_trees_.cbegin(); }

        AnnotatedTreeIterator cend() const { return annotated_trees_.cend(); }

        size_t size() const { return annotated_trees_.size(); }
    };
}