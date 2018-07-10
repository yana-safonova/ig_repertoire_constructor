#pragma once

#include "clonal_graph.hpp"
#include "aa_shm_map.hpp"

namespace antevolo {
    class ClonalTreeRefiner {
        typedef std::pair<size_t, size_t> ClonalEdge;

        const ClonalGraph& clonal_graph_;

        AaSHMMap aa_shm_map_;
        EvolutionaryTree refined_tree_;
        /*
        void FillAaSHMMap();

        ClonalTreeRefiner::ClonalEdge GetMaximalConflictingEdge();

        void ResolveConflictingGroup(ClonalEdge selected_edge);

        bool ConflictsExist();
        */

    public:
        ClonalTreeRefiner(const ClonalGraph& clonal_graph,
                          CloneSetWithFakesPtr clone_set_ptr) : clonal_graph_(clonal_graph),
                                                                aa_shm_map_(),
                                                                refined_tree_(clone_set_ptr) { }
        /*
        EvolutionaryTree RefineClonalTree() {
            while(ConflictsExist()) {
                ClonalEdge edge = GetMaximalConflictingEdge();
                ResolveConflictingGroup(edge);
            }
            return refined_tree_;
        }
        */
    };
}