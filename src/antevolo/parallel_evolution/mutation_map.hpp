#pragma once

#include "clonal_graph.hpp"
#include "../shm_counting/tree_shm_calculator.hpp"
#include "../clone_set_with_fakes.hpp"

namespace antevolo {
    class MutationMap {
        const CloneSetWithFakes &clone_set_;
        //const ParallelEvolutionStats &stats_;
        const ClonalGraph &clonal_graph_;

        std::map<std::pair<size_t, size_t>, std::vector<TreeSHM> > edge_shms_map_;
        std::map<TreeSHM, size_t> unique_shm_mult_;
        std::map<TreeSHM, std::vector<std::pair<size_t, size_t> > > shm_edge_map_;

        void FillMutationFromEdge(size_t src, size_t dst);

        void FillMutationsFromGraph();

        void FillSHMEdgeMap();

        void ComputeUniqueMutationMult();

        //void FindUntrustedMutations();

        bool SHMIsHotspot(TreeSHM shm, std::pair<size_t, size_t> edge) const;

    public:
        MutationMap(const CloneSetWithFakes &clone_set,
                    const ClonalGraph &clonal_graph) : clone_set_(clone_set),
                                                       clonal_graph_(clonal_graph) {
            FillMutationsFromGraph();
            FillSHMEdgeMap();
            ComputeUniqueMutationMult();
        }

        size_t NumSynonymousNonTrivialSHMs() const;

        size_t NumNontrivialSHMs() const;

        size_t MaxMultiplicity() const;

        size_t NumNonTrivialHotSpotSHMs() const;

        size_t NumUniqueSHMs() const { return unique_shm_mult_.size(); }


        bool SHMIsHotSpot(TreeSHM shm) const;

        bool SHMIsNonTrivial(TreeSHM shm) const;

        size_t GetSHMMultiplicity(TreeSHM shm) const;


        typedef std::map<TreeSHM, size_t>::const_iterator SHMMultiplicityConstIterator;

        SHMMultiplicityConstIterator shm_mult_cbegin() const { return unique_shm_mult_.cbegin(); }

        SHMMultiplicityConstIterator shm_mult_cend() const { return unique_shm_mult_.cend(); }


        typedef std::map<std::pair<size_t, size_t>, std::vector<TreeSHM> >::const_iterator EdgeSHMsMapConstIterator;

        EdgeSHMsMapConstIterator edge_shm_cbegin() const { return edge_shms_map_.cbegin(); }

        EdgeSHMsMapConstIterator edge_shm_cend() const { return edge_shms_map_.cend(); }

        std::vector<TreeSHM> GetSHMsByEdge(size_t src, size_t dst) const;

        std::vector<std::pair<size_t, size_t> > GetEdgesBySHM(TreeSHM shm) const;


        const CloneSetWithFakes& CloneSet() const { return clone_set_; }

        const ClonalGraph& GetClonalGraph() const { return clonal_graph_; }

        seqan::CharString VGeneName() const;


        std::set<TreeSHM> GetUniqueOutgoingSHMs(size_t src) const;
    };
}