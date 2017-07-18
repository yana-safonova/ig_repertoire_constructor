#pragma once

#include "tree_based_shm.hpp"
#include "../evolutionary_graph_utils/evolutionary_tree.hpp"

namespace antevolo {
    class TreeSHMMap {
        //const CloneSetWithFakes &clone_set_;
        //const EvolutionaryTree &tree_;

        std::map<TreeSHM, size_t> shm_mult_map_;
        std::map<TreeSHM, std::vector<std::pair<size_t, size_t> > > shm_clone_ids_;

    public:
        TreeSHMMap(/*const CloneSetWithFakes &clone_set,
                   const EvolutionaryTree &tree*/) { } //: clone_set_(clone_set), tree_(tree) { }

        void AddSHM(TreeSHM shm, size_t src, size_t dst_id);

        typedef std::map<TreeSHM, size_t>::const_iterator TreeSHMConstIterator;

        TreeSHMConstIterator cbegin() const { return shm_mult_map_.cbegin(); }

        TreeSHMConstIterator cend() const { return shm_mult_map_.cend(); }

        typedef std::map<TreeSHM, std::vector<std::pair<size_t, size_t> > >::const_iterator TreeSHMCloneConstIterator;

        TreeSHMCloneConstIterator c_shm_clone_begin() const { return shm_clone_ids_.cbegin(); }

        TreeSHMCloneConstIterator c_shm_clone_end() const { return shm_clone_ids_.cend(); }

        size_t size() const { return shm_mult_map_.size(); }

        size_t NumSynonymousSHMs() const;

        size_t NumCDRSHMs() const;

        size_t MaxMultiplicity() const;

        size_t NumSHMsInRegion(annotation_utils::StructuralRegion region) const;

        size_t NumSynSHMsInRegion(annotation_utils::StructuralRegion region) const;
    };
}