#include "sparse_graph.hpp"

// vertex set should be sorted
std::shared_ptr<SparseGraph> SparseGraph::GetSubgraph(size_t subgraph_id, const set<size_t> &vertex_set) {
    vector<GraphEdge> subgraph_edges;
    component_map_.AddComponentInMap(subgraph_id, vertex_set);
    for(auto it = vertex_set.begin(); it != vertex_set.end(); it++) {
        size_t vertex1 = *it;
        for(size_t i = RowIndex()[vertex1]; i < RowIndex()[vertex1 + 1]; i++) {
            size_t vertex2 = Col()[i];
            size_t weight = Dist()[i];
            if(vertex_set.find(vertex2) != vertex_set.end()) {
                size_t new_vertex1 = component_map_.GetNewVertexByOldVertex(vertex1);
                size_t new_vertex2 = component_map_.GetNewVertexByOldVertex(vertex2);
                subgraph_edges.push_back(GraphEdge(new_vertex1, new_vertex2, weight));
            }
        }
    }
    INFO("Subgraph contains " << vertex_set.size() << " vertices and " << subgraph_edges.size() << " edges");
    return std::shared_ptr<SparseGraph>(new SparseGraph(vertex_set.size(), subgraph_edges));
}

ostream& operator<<(ostream &out, const SparseGraph &graph) {
    out << "Direct matrix" << endl;
    out << *(graph.DirectMatrix()) << endl;
    out << "-------------" << endl;
    out << "Transposed matrix" << endl;
    out << *(graph.TransposedMatrix());
    return out;
}