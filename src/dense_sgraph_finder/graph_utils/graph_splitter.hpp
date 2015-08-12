#include "sparse_graph.hpp"

class ConnectedComponentGraphSplitter {
    SparseGraphPtr graph_ptr_;
    map<size_t, set<size_t>> neighbourhood_map_;
    set<size_t> visited_vertices_;

    void InitializeMap();

    size_t GetStartVertex();

    SparseGraphPtr GetConnectedComponentByVertex(size_t start_vertex);

public:
    ConnectedComponentGraphSplitter(SparseGraphPtr graph_ptr) : graph_ptr_(graph_ptr) { }

    vector<SparseGraphPtr> Split();
};