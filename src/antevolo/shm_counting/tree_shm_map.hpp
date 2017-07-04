#pragma once

#include "tree_based_shm.hpp"
#include "../evolutionary_graph_utils/evolutionary_tree.hpp"

namespace antevolo {
    class TreeSHMMap {
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;
        const EvolutionaryTree &tree_;

        std::map<TreeSHM, size_t> shm_mult_map_;

    public:
        TreeSHMMap(const annotation_utils::CDRAnnotatedCloneSet &clone_set,
                   const EvolutionaryTree &tree) : clone_set_(clone_set), tree_(tree) { }

        void AddSHM(TreeSHM shm);

        typedef std::map<TreeSHM, size_t>::const_iterator TreeSHMIterator;

        TreeSHMIterator cbegin() const { return shm_mult_map_.cbegin(); }

        TreeSHMIterator cend() const { return shm_mult_map_.cend(); }

        size_t size() const { return shm_mult_map_.size(); }
    };
}