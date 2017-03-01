#pragma once

#include "evolutionary_graph_utils/evolutionary_tree.hpp"

namespace antevolo {
    class EvolutionaryTreeStorage {
//        const annotation_utils::CDRAnnotatedCloneSet* clone_set_ptr_;

        std::vector<EvolutionaryTree> trees_;

    public:
//        EvolutionaryTreeStorage(const annotation_utils::CDRAnnotatedCloneSet& clone_set) :
//                clone_set_ptr_(&clone_set) { }

        void Add(EvolutionaryTree tree);

        void AppendArchive(const EvolutionaryTreeStorage& tree_storage);

        typedef std::vector<EvolutionaryTree>::iterator TreeStorageIterator;

        typedef std::vector<EvolutionaryTree>::const_iterator TreeStorageConstIterator;

        TreeStorageIterator begin() { return trees_.begin(); }

        TreeStorageIterator end() { return trees_.end(); }

        TreeStorageConstIterator cbegin() const { return trees_.cbegin(); }

        TreeStorageConstIterator cend() const { return trees_.cend(); }

        size_t size() const { return trees_.size(); }

        const EvolutionaryTree& GetTreeByIndex(size_t index) const;

        const EvolutionaryTree& operator[](size_t index) const;
    };
}