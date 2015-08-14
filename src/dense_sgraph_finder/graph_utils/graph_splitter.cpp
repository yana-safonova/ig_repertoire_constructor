#include "graph_splitter.hpp"

void ConnectedComponentGraphSplitter::InitializeMap() {
    for(size_t i = 0; i < graph_ptr_->N(); i++) {
        for(size_t j = graph_ptr_->RowIndex()[i]; j < graph_ptr_->RowIndex()[i + 1]; j++) {
            size_t col = graph_ptr_->Col()[j];
            if (neighbourhood_map_.find(i) == neighbourhood_map_.end())
                neighbourhood_map_[i] = set<size_t>();
            neighbourhood_map_[i].insert(col);
            if (neighbourhood_map_.find(col) == neighbourhood_map_.end())
                neighbourhood_map_[col] = set<size_t>();
            neighbourhood_map_[col].insert(i);
        }
    }
}

size_t ConnectedComponentGraphSplitter::GetStartVertex() {
    for(size_t i = 0; i < graph_ptr_->N(); i++)
        if(visited_vertices_.find(i) == visited_vertices_.end())
            return i;
    return size_t(-1);
}


SparseGraphPtr ConnectedComponentGraphSplitter::GetConnectedComponentByVertex(size_t component_id, size_t start_vertex) {
    queue<size_t> vertex_queue;
    vertex_queue.push(start_vertex);
    set<size_t> connected_component;
    while(!vertex_queue.empty()) {
        size_t cur_vertex = vertex_queue.front();
        vertex_queue.pop();
        connected_component.insert(cur_vertex);
        visited_vertices_.insert(cur_vertex);
        for(size_t i = graph_ptr_->RowIndex()[cur_vertex]; i < graph_ptr_->RowIndex()[cur_vertex + 1]; i++) {
            size_t cur_neigh = graph_ptr_->Col()[i];
            if(visited_vertices_.find(cur_neigh) == visited_vertices_.end()) {
                vertex_queue.push(cur_neigh);
                visited_vertices_.insert(cur_neigh);
            }
        }
        for(size_t i = graph_ptr_->RowIndexT()[cur_vertex]; i < graph_ptr_->RowIndexT()[cur_vertex + 1]; i++) {
            size_t cur_neigh = graph_ptr_->ColT()[i];
            if(visited_vertices_.find(cur_neigh) == visited_vertices_.end()) {
                vertex_queue.push(cur_neigh);
                visited_vertices_.insert(cur_neigh);
            }
        }
    }
    return graph_ptr_->GetSubgraph(component_id, connected_component);
}

vector<SparseGraphPtr> ConnectedComponentGraphSplitter::Split() {
    InitializeMap();
    INFO("Connected component splitter was initialized");
    vector<SparseGraphPtr> connected_components;
    size_t start_vertex = GetStartVertex();
    while(start_vertex != size_t(-1)) {
        size_t component_id = connected_components.size();
        SparseGraphPtr cur_connected_component = GetConnectedComponentByVertex(component_id, start_vertex);
        connected_components.push_back(cur_connected_component);
        start_vertex = GetStartVertex();
    }
    INFO("Graph was splitted into " << connected_components.size() << " connected component(s)");
    return connected_components;
}