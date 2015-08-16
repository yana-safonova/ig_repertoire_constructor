#include "graph_splitter.hpp"

void ConnectedComponentGraphSplitter::InitializeInnerVertices() {
    for(size_t i = 0; i < graph_ptr_->N(); i++)
        visited_vertices_.push_back(false);
}

size_t ConnectedComponentGraphSplitter::GetStartVertex() {
    for(size_t i = next_start_vertex_; i < graph_ptr_->N(); i++)
        if(!visited_vertices_[i]) {
            next_start_vertex_ = i + 1;
            return i;
        }
    return size_t(-1);
}


SparseGraphPtr ConnectedComponentGraphSplitter::GetConnectedComponentByVertex(size_t component_id, size_t start_vertex) {
    queue<size_t> vertex_queue;
    vertex_queue.push(start_vertex);
    set<size_t> connected_component;
    // check if vertex is isolated
    if(graph_ptr_->VertexIsIsolated(start_vertex)) {
        connected_component.insert(start_vertex);
        return graph_ptr_->GetSubgraph(component_id, connected_component);
    }
    // otherwise search for connected component
    while(!vertex_queue.empty()) {
        size_t cur_vertex = vertex_queue.front();
        vertex_queue.pop();
        connected_component.insert(cur_vertex);
        visited_vertices_[cur_vertex] = true;
        for(size_t i = graph_ptr_->RowIndex()[cur_vertex]; i < graph_ptr_->RowIndex()[cur_vertex + 1]; i++) {
            size_t cur_neigh = graph_ptr_->Col()[i];
            if(!visited_vertices_[cur_neigh]) {
                vertex_queue.push(cur_neigh);
                visited_vertices_[cur_neigh] = true;
            }
        }
        for(size_t i = graph_ptr_->RowIndexT()[cur_vertex]; i < graph_ptr_->RowIndexT()[cur_vertex + 1]; i++) {
            size_t cur_neigh = graph_ptr_->ColT()[i];
            if(!visited_vertices_[cur_neigh]) {
                vertex_queue.push(cur_neigh);
                visited_vertices_[cur_neigh] = true;
            }
        }
    }
    return graph_ptr_->GetSubgraph(component_id, connected_component);
}

void ConnectedComponentGraphSplitter::PrintConnectedComponentsStats(const vector<SparseGraphPtr> &connected_components) {
    size_t max_vertex_size = 0;
    size_t max_edge_size = 0;
    size_t num_small_components = 0;
    size_t num_singletons = 0;
    // todo: move it to the config
    size_t min_graph_size = 4;
    for(auto it = connected_components.begin(); it != connected_components.end(); it++) {
        if((*it)->N() == 1)
            num_singletons++;
        if((*it)->N() <= min_graph_size)
            num_small_components++;
        if((*it)->N() > max_vertex_size) {
            max_vertex_size = (*it)->N();
            max_edge_size = (*it)->NZ();
        }
        else if((*it)->N() == max_vertex_size and (*it)->NZ() > max_edge_size)
            max_edge_size = (*it)->NZ();
    }
    INFO("Largest component contains " << max_vertex_size << " vertices & " << max_edge_size << " edges");
    float singleton_perc = (float)num_singletons / float(connected_components.size());
    INFO("# singleton components: " << num_singletons << " (" << singleton_perc << "%)");
    float small_comp_perc = float(num_small_components) / float(connected_components.size());
    INFO("# small components (# vertices <= " << min_graph_size << "): " <<
                 num_small_components << " (" << small_comp_perc << "%)");
}

vector<SparseGraphPtr> ConnectedComponentGraphSplitter::Split() {
    InitializeInnerVertices();
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
    PrintConnectedComponentsStats(connected_components);
    return connected_components;
}