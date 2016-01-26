#include <verify.hpp>
#include "sparse_graph.hpp"


size_t SparseGraph::get_edge_at_index(size_t vertex, size_t idx) const {
    return idx < RowIndexT()[vertex + 1] - RowIndexT()[vertex]
           ? ColT()[RowIndexT()[vertex] + idx]
           : Col()[RowIndex()[vertex] + idx - (RowIndexT()[vertex + 1] - RowIndexT()[vertex])];
}

bool SparseGraph::HasEdge(size_t from, size_t to) const {
    // for bidirectional graphs only (we support only them currently)
    size_t from_deg = Degree(from);
    size_t to_deg = Degree(to);
    if (to_deg < from_deg) return HasEdge(to, from);

    if (from_deg == 0) return false;

    size_t left = 0;
    size_t right = from_deg - 1;
    while (left < right) {
        size_t mid = (left + right) / 2;
        size_t v = get_edge_at_index(from, mid);
        if (v < to) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    return to == get_edge_at_index(from, right);
}

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
    //INFO("Subgraph contains " << vertex_set.size() << " vertices and " << subgraph_edges.size() << " edges");
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

SparseGraph::EdgesIterator SparseGraph::EdgesIterator::operator++() {
    current_++;
    return *this;
}

SparseGraph::EdgesIterator SparseGraph::EdgesIterator::operator++(int) {
    const SparseGraph::EdgesIterator& itr = SparseGraph::EdgesIterator(*this);
    current_++;
    return itr;
}

bool SparseGraph::EdgesIterator::operator==(SparseGraph::EdgesIterator other) const {
    VERIFY_MSG(&vertex_ == &other.vertex_, "Comparing iterators over different vertices edges.");
    return current_ == other.current_;
}

bool SparseGraph::EdgesIterator::operator!=(SparseGraph::EdgesIterator other) const {
    VERIFY_MSG(&vertex_ == &other.vertex_, "Comparing iterators over different vertices edges.");
    return current_ != other.current_;
}

size_t SparseGraph::EdgesIterator::operator*() const {
    return graph_.get_edge_at_index(vertex_.GetIndex(), current_);
}
