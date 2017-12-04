#pragma once

#include "tree_shm_map.hpp"
#include "../evolutionary_graph_utils/evolutionary_tree.hpp"
#include "../clone_set_with_fakes.hpp"

namespace antevolo {
    class TreeSHMCalculator {
        const CloneSetWithFakes &clone_set_;

    public:
        TreeSHMCalculator(const CloneSetWithFakes &clone_set) : clone_set_(clone_set) { }

        TreeSHM ComputeTreeSHMByUsualSHM(annotation_utils::SHM shm, size_t src, size_t dst) const;
    };

    class UniqueSHMCalculator {
        const CloneSetWithFakes &clone_set_;
        const EvolutionaryTree &tree_;

        TreeSHMCalculator tree_shm_calc_;
        TreeSHMMap shm_map_;

        void AddAddedSHMs(std::vector<annotation_utils::SHM> &shms, size_t src, size_t dst);

        void AddNestedSHMs(size_t src, size_t dst);

        void AddCDR3SHMs(size_t src, size_t dst);

        bool FilterEdge(size_t src, size_t dst) const;

        TreeSHM ComputeSHMInCDR3(size_t src, size_t dst, size_t cdr3_pos) const;

    public:
        UniqueSHMCalculator(const CloneSetWithFakes &clone_set,
                            const EvolutionaryTree &tree) : clone_set_(clone_set),
                                                            tree_(tree),
                                                            tree_shm_calc_(clone_set_),
                                                            shm_map_(){ //(clone_set_, tree) {
        }

        void AddSHMsFromEdge(const EvolutionaryEdgePtr edge_ptr);

        void AddSHMsFromEdge(size_t src_id, size_t dst_id);

        TreeSHMMap GetSHMMap() { return shm_map_; }
    };
}
