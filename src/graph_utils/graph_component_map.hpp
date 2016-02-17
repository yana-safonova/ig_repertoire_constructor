#pragma once

#include "include_me.hpp"

class GraphComponentMap {
    map<size_t, size_t> old_vertex_to_subgraph_;
    map<size_t, size_t> old_vertex_to_new_vertex_;
    map<size_t, map<size_t, size_t>> component_map_;
    vector<size_t> subgraph_ids_;
    vector<size_t> old_vertices_list_;

    void InitializeSubgraphIds();

    void InitializeOldVerticesList();

public:
    GraphComponentMap() { }

    void AddComponentInMap(size_t subgraph_id, const set<size_t> &old_vertices_set);

    size_t GetSubgraphIdByOldVertex(size_t old_vertex);

    size_t GetNewVertexByOldVertex(size_t old_vertex);

    size_t GetOldVertexByNewVertex(size_t subgraph_id, size_t new_vertex);

    const vector<size_t>& SubgraphIds();

    const vector<size_t>& OldVerticesList();
};

ostream& operator<<(ostream &out, GraphComponentMap& component_map);

typedef shared_ptr<GraphComponentMap> GraphComponentMapPtr;