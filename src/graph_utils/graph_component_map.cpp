#include <verify.hpp>
#include "graph_component_map.hpp"

void GraphComponentMap::AddComponentInMap(size_t subgraph_id, const std::set<size_t> &old_vertices_set) {
    VERIFY(component_map_.find(subgraph_id) == component_map_.end());
    std::map<size_t, size_t> cur_component_map;
    size_t new_vertex_id = 0;
    for(auto it = old_vertices_set.begin(); it != old_vertices_set.end(); it++) {
        cur_component_map[new_vertex_id] = *it;
        old_vertex_to_new_vertex_[*it] = new_vertex_id;
        old_vertex_to_subgraph_[*it] = subgraph_id;
        new_vertex_id++;
    }
    component_map_[subgraph_id] = cur_component_map;
}

size_t GraphComponentMap::GetSubgraphIdByOldVertex(size_t old_vertex) {
    VERIFY(old_vertex_to_subgraph_.find(old_vertex) != old_vertex_to_subgraph_.end());
    return old_vertex_to_subgraph_[old_vertex];
}

size_t GraphComponentMap::GetNewVertexByOldVertex(size_t old_vertex) {
    VERIFY(old_vertex_to_new_vertex_.find(old_vertex) != old_vertex_to_new_vertex_.end());
    return old_vertex_to_new_vertex_[old_vertex];
}

void GraphComponentMap::InitializeSubgraphIds() {
    for(auto it = component_map_.begin(); it != component_map_.end(); it++)
        subgraph_ids_.push_back(it->first);
}

size_t GraphComponentMap::GetOldVertexByNewVertex(size_t subgraph_id, size_t new_vertex) {
    VERIFY(component_map_.find(subgraph_id) != component_map_.end());
    std::map<size_t, size_t> &cur_component_map = component_map_[subgraph_id];
    VERIFY(cur_component_map.find(new_vertex) != cur_component_map.end());
    return cur_component_map[new_vertex];
}

const std::vector<size_t>& GraphComponentMap::SubgraphIds() {
    if(subgraph_ids_.size() == 0)
        InitializeSubgraphIds();
    return subgraph_ids_;
}

void GraphComponentMap::InitializeOldVerticesList() {
    for(auto it = old_vertex_to_subgraph_.begin(); it != old_vertex_to_subgraph_.end(); it++) {
        VERIFY(old_vertex_to_new_vertex_.find(it->first) != old_vertex_to_new_vertex_.end());
        old_vertices_list_.push_back(it->first);
    }
}

const std::vector<size_t>& GraphComponentMap::OldVerticesList() {
    if(old_vertices_list_.size() == 0)
        InitializeOldVerticesList();
    return old_vertices_list_;
}

std::ostream& operator<<(std::ostream &out, GraphComponentMap& component_map) {
    out << "Old vertex -> subgraph id, new vertex" << std::endl;
    for(auto it = component_map.OldVerticesList().cbegin(); it != component_map.OldVerticesList().cend(); it++) {
        out << *it << " -> " << component_map.GetSubgraphIdByOldVertex(*it) << ", " <<
                component_map.GetNewVertexByOldVertex(*it) << std::endl;
    }
    return out;
}