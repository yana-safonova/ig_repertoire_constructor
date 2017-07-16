#pragma once

#include "parallel_rhomb.hpp"
#include "../clone_set_with_fakes.hpp"

namespace antevolo {
    struct ParallelEvolutionStats {
        size_t num_parallel_rhombs;
        size_t total_num_parallel_shms;
        size_t num_ambiguous_rhombs;
        std::vector<size_t> parallel_shms_dist;
        std::vector<ParallelRhomb> rhombs;

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
            rhombs.push_back(rhomb);
        }

        std::vector<size_t>::iterator begin() { return parallel_shms_dist.begin(); }

        std::vector<size_t>::iterator end() { return parallel_shms_dist.end(); }

        void AddParallelSHM(size_t num_parallel_shms) {
            parallel_shms_dist.push_back(num_parallel_shms);
            total_num_parallel_shms += num_parallel_shms;
        }

        bool Empty() const { return rhombs.size() == 0; }
    };

    class ParallelEvolutionFinder {
        const EvolutionaryTree &tree_;
        const CloneSetWithFakes &clone_set_;

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
        ParallelEvolutionFinder(const EvolutionaryTree &tree) : tree_(tree),
                                                                clone_set_(*tree.GetCloneSetPtr()) { }

        ParallelEvolutionStats ComputeParallelSHMs();
    };
}
