#pragma once

#include "sparse_graph.hpp"

class ConnectedComponentGraphSplitter {
    // input parameters
    SparseGraphPtr graph_ptr_;

    // auxiliary parameters
    vector<bool> visited_vertices_;
    size_t next_start_vertex_;

    void InitializeInnerVertices();

    size_t GetStartVertex();

    SparseGraphPtr GetConnectedComponentByVertex(size_t component_id, size_t start_vertex);

public:
    ConnectedComponentGraphSplitter(SparseGraphPtr graph_ptr) :
            graph_ptr_(graph_ptr),
            next_start_vertex_(0) { }

    vector<SparseGraphPtr> Split();
};