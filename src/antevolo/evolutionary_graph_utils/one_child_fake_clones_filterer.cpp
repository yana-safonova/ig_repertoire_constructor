#include "one_child_fake_clones_filterer.hpp"

namespace antevolo {
    EvolutionaryTree OneChildFakeClonesFilterer::FilterOneChildFakes(const EvolutionaryTree& connected_tree) const {
//        tree have to be connected!!
        EvolutionaryTree filtered_tree(connected_tree.GetCloneSetPtr());
        auto edge_constructor = GetEdgeConstructor();
        const auto& clone_set = connected_tree.GetCloneSet();

        size_t root = connected_tree.GetRoots()[0];
        VERIFY(connected_tree.GetRoots().size() == 1);
        std::queue<size_t> vertex_queue;
        vertex_queue.push(root);
        filtered_tree.SetTreeIndices(connected_tree.GetVJClassIndex(),
                                     connected_tree.GetConnectedComponentIndex(),
                                     connected_tree.GetTreeIndex());
        while(!vertex_queue.empty()) {
            size_t cur_vertex = vertex_queue.front();
            vertex_queue.pop();
            if(!connected_tree.IsRoot(cur_vertex) && !connected_tree.IsFakeToFilter(cur_vertex)) {
                EvolutionaryEdgePtr edge = connected_tree.GetParentEdge(cur_vertex);
                size_t current_parent = edge->SrcNum();
                while (current_parent != root &&
                       connected_tree.IsFakeToFilter(current_parent)) {
                    current_parent = connected_tree.GetParentEdge(current_parent)->SrcNum();
                }
                if (!connected_tree.IsFakeToFilter(current_parent)) {
                    edge = edge_constructor->ConstructEdge(clone_set[current_parent],
                                                           clone_set[cur_vertex],
                                                           current_parent,
                                                           cur_vertex);
                    VERIFY_MSG(edge->IsDirected() ||
                               edge->IsUndirected() ||
                               edge->IsReverseDirected() ||
                               edge->IsDoubleMutated() ||
                               edge->IsIntersected(),
                               "Edge from " << edge->DstNum() << " -> " << edge->SrcNum() <<
                                            " is not directed, undirected, reverse directed, double mutated "
                                            << "or even intersected");
                    filtered_tree.AddEdge(edge->DstNum(), edge);
                }
            }
            if(!connected_tree.IsLeaf(cur_vertex) ) {
                const auto& outgoing_edges = connected_tree.OutgoingEdges(cur_vertex);
                for(auto it = outgoing_edges.begin(); it != outgoing_edges.end(); it++) {
                    EvolutionaryEdgePtr edge = *it;
                    vertex_queue.push(edge->DstNum());
                }
            }
        }
        return filtered_tree;
    }
}