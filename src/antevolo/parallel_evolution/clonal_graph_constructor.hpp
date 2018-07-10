#pragma once

#include "clonal_graph.hpp"

namespace antevolo {
    class ClonalGraphConstructor {
        const EvolutionaryTree &tree_;
        const CloneSetWithFakes &clone_set_;

        ClonalGraph clonal_graph_;

        void AddEdgesFromTree();

        void AddDirectedEdges();

        //void InitializeClonalGraph();

        bool SkipVertex(size_t vertex) const;

        bool VerticesConnectedByDirectedPath(size_t v_src, size_t v_dst) const;

        bool VerticesCanBeConnected(size_t v_src, size_t v_dst) const;

        bool VerticesBelongToSameTree(size_t v_src, size_t v_dst) const;

        bool ParentEdgeIsUndirected(size_t node) const;

    public:
        ClonalGraphConstructor(const EvolutionaryTree &tree) : tree_(tree),
                                                               clone_set_(tree.GetCloneSet()),
                                                               clonal_graph_(/*tree*/) { }
        ClonalGraph ComputeClonalGraph() {
            AddEdgesFromTree();
            AddDirectedEdges();
            clonal_graph_.RemoveExtraEdges();
            return clonal_graph_;
        }
    };
}