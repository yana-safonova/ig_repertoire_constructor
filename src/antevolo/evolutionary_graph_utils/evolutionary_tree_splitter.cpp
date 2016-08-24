#include "evolutionary_tree_splitter.hpp"

namespace antevolo {
    EvolutionaryTree ConnectedTreeSplitter::GetTreeByRoot(const EvolutionaryTree &tree,
                                                                           size_t root_id) {
        std::queue<size_t> vertex_queue;
        vertex_queue.push(root_id);
        EvolutionaryTree connected_tree;
        while(!vertex_queue.empty()) {
            size_t cur_vertex = vertex_queue.front();
            vertex_queue.pop();
            if(!tree.IsRoot(cur_vertex)) {
                EvolutionaryEdge edge = tree.GetParentEdge(cur_vertex);
                // todo: refactor it!
                if(edge.IsDirected()) {
                    connected_tree.AddDirected(edge.dst_clone_num, edge);
                }
                else if(edge.IsUndirected()) {
                    connected_tree.AddUndirected(edge.dst_clone_num, edge);
                }
                else {
                    VERIFY_MSG(false, "Edge from " << edge.dst_clone_num << " -> " << edge.src_clone_num <<
                            " is not directed or undirected");
                }
            }
            if(!tree.IsLeaf(cur_vertex)) {
                auto outgoing_edges = tree.OutgoingEdges(cur_vertex);
                for(auto it = outgoing_edges.begin(); it != outgoing_edges.end(); it++) {
                    EvolutionaryEdge edge = *it;
                    vertex_queue.push(edge.dst_clone_num);
                }
            }
        }
        return connected_tree;
    }

    std::vector<EvolutionaryTree> ConnectedTreeSplitter::Split(const EvolutionaryTree &tree) {
        std::vector<EvolutionaryTree> connected_trees;
        if(tree.GetRootNumber() == 1) {
            connected_trees.push_back(tree);
            return connected_trees;
        }
        auto roots = tree.GetRoots();
        for(auto it = roots.begin(); it != roots.end(); it++) {
            EvolutionaryTree connected_tree = GetTreeByRoot(tree, *it);
            connected_trees.push_back(connected_tree);
        }
        VERIFY_MSG(connected_trees.size() == tree.GetRootNumber(), "ERROR: number of connected trees (" <<
                connected_trees.size() << " does not match with number of roots (" << tree.GetRootNumber() <<
                ") in the original forest");
        return connected_trees;
    }
}