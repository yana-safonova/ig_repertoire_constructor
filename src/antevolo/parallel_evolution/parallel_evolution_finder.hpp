#pragma once

#include "parallel_rhomb.hpp"

namespace antevolo {
    struct ParallelEvolutionStats {
        size_t num_parallel_rhombs;
        size_t total_num_parallel_shms;
        size_t num_ambiguous_rhombs;
        std::vector<size_t> parallel_shms_dist;

        ParallelEvolutionStats() : num_parallel_rhombs(0),
                                   total_num_parallel_shms(0),
                                   num_ambiguous_rhombs(0) { }

        void ConcatenateStats(ParallelEvolutionStats stats) {
            num_parallel_rhombs += stats.num_parallel_rhombs;
            num_ambiguous_rhombs += stats.num_ambiguous_rhombs;
            for(auto it = stats.begin(); it != stats.end(); it++)
                AddParallelSHM(*it);
        }

        void Add(const ParallelRhomb &rhomb) {
            num_parallel_rhombs++;
            if(rhomb.Ambiguous()) {
                num_ambiguous_rhombs++;
            }
            AddParallelSHM(rhomb.MinimalNumberParallelSHMs());
        }

        std::vector<size_t>::iterator begin() { return parallel_shms_dist.begin(); }

        std::vector<size_t>::iterator end() { return parallel_shms_dist.end(); }

        void AddParallelSHM(size_t num_parallel_shms) {
            parallel_shms_dist.push_back(num_parallel_shms);
            total_num_parallel_shms += num_parallel_shms;
        }
    };

    class ParallelEvolutionFinder {
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;
        const EvolutionaryTree &tree_;

        std::map<size_t, std::set<size_t>> added_directed_edges_map_;

        bool VerticesConnectedByDirectedPath(size_t v_src, size_t v_dst) const;

        bool VerticesBelongToSameTree(size_t v_src, size_t v_dst) const;

        bool ParentEdgeIsUndirected(size_t node) const;

        bool PairIsGood(size_t v_src, size_t v_dst) const;

        void FillAddedDirectedEdges();

        void RemoveUndirectedRhombs();

        void RemoveIncomingNestedRhombs();

        void RemoveNestedOutgoingRhombs();

        void PrintAddedNestedEdges();

        EvolutionaryEdgePtr GetAddedEvolutionaryTree(size_t v_src, size_t v_dst);

        ParallelRhomb ComputeRhombByAddedEdge(size_t v_src, size_t v_dst);

        bool RhombIsGood(ParallelRhomb rhomb);

        std::vector<ParallelRhomb> ComputeRhombs();

        ParallelEvolutionStats ComputeStats(const std::vector<ParallelRhomb>& rhombs);

        //bool VerticesPresentTrueParallelEvolution(size_t src_id, size_t dst_id);

    public:
        ParallelEvolutionFinder(const annotation_utils::CDRAnnotatedCloneSet& clone_set,
                                const EvolutionaryTree &tree) : clone_set_(clone_set),
                                                                tree_(tree) { }

        ParallelEvolutionStats ComputeParallelSHMs();
    };
}
