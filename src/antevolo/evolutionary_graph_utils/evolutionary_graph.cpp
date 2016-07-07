#include <logger/logger.hpp>
#include "evolutionary_graph.hpp"

namespace antevolo {
    size_t EvolutionaryGraph::GetVertexIndex(const annotation_utils::AnnotatedClone &clone) {
        if(vertex_index_map_.find(clone.Read().name) != vertex_index_map_.end())
            return vertex_index_map_[clone.Read().name];
        size_t map_size = vertex_index_map_.size();
        vertex_index_map_[clone.Read().name] =  map_size;
        return map_size;
    }

    void EvolutionaryGraph::UpdateVertexMap(size_t src_index, size_t dst_index) {
        if(vertex_neighbours_map_.find(src_index) == vertex_neighbours_map_.end())
            vertex_neighbours_map_[src_index] = std::vector<size_t>();
        vertex_neighbours_map_[src_index].push_back(dst_index);
    }

    void EvolutionaryGraph::UpdateEdgeMap(size_t vertex_index, size_t edge_index) {
        if(vertex_edges_map_.find(vertex_index) == vertex_edges_map_.end())
            vertex_edges_map_[vertex_index] = std::vector<size_t>();
        vertex_edges_map_[vertex_index].push_back(edge_index);
    }

    void EvolutionaryGraph::AddEdge(EvolutionaryEdge edge) {
        if(edge.Empty())
            return;
        edges_.push_back(edge);
        size_t edge_index = edges_.size() - 1;
        if(edge_type_counts_.find(edge.edge_type) == edge_type_counts_.end())
            edge_type_counts_[edge.edge_type] = 0;
        edge_type_counts_[edge.edge_type]++;
        // add new vertex
        size_t src_clone_index = GetVertexIndex(*edge.src_clone);
        size_t dst_clone_index = GetVertexIndex(*edge.dst_clone);
        UpdateEdgeMap(src_clone_index, edge_index);
        UpdateVertexMap(src_clone_index, dst_clone_index);
        if(!edge.IsDirected()) {
            UpdateEdgeMap(dst_clone_index, edge_index);
            UpdateVertexMap(dst_clone_index, src_clone_index);
        }
    }

    size_t EvolutionaryGraph::NumEdgesOfType(EvolutionaryEdgeType edge_type) const {
        if(edge_type_counts_.find(edge_type) == edge_type_counts_.end())
            return 0;
        return edge_type_counts_.at(edge_type);
    }

    EvolutionaryTree EvolutionaryGraph::GetEvolutionaryTree() {
        std::priority_queue<EvolutionaryEdge, std::vector<EvolutionaryEdge>,
                WeightEvolutionaryEdgeComparator> edge_queue;
        for(auto it = edges_.begin(); it != edges_.end(); it++)
            edge_queue.push(*it);
        EvolutionaryTree tree;
        std::set<std::string> clone_set;
        while(!edge_queue.empty()) {
            EvolutionaryEdge edge = edge_queue.top();
            if(clone_set.find(edge.dst_clone->Read().name) == clone_set.end() or
                    clone_set.find(edge.src_clone->Read().name) == clone_set.end()) {
                tree.Add(edge);
                clone_set.insert(edge.dst_clone->Read().name);
                clone_set.insert(edge.src_clone->Read().name);
            }
            edge_queue.pop();
        }
        TRACE("Evolutionary tree contains " << tree.NumEdges() << " edges");
        return tree;
    }
}