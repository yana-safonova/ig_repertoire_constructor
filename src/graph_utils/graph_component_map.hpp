#pragma once

#include <map>
#include <vector>
#include <set>
#include <memory>

class GraphComponentMap {
    std::map<size_t, size_t> old_vertex_to_subgraph_;
    std::map<size_t, size_t> old_vertex_to_new_vertex_;
    std::map<size_t, std::map<size_t, size_t>> component_map_;
    std::vector<size_t> subgraph_ids_;
    std::vector<size_t> old_vertices_list_;

    void InitializeSubgraphIds();

    void InitializeOldVerticesList();

public:
    GraphComponentMap() { }

    void AddComponentInMap(size_t subgraph_id, const std::set<size_t> &old_vertices_set);

    size_t GetSubgraphIdByOldVertex(size_t old_vertex);

    size_t GetNewVertexByOldVertex(size_t old_vertex);

    size_t GetOldVertexByNewVertex(size_t subgraph_id, size_t new_vertex);

    const std::vector<size_t>& SubgraphIds();

    const std::vector<size_t>& OldVerticesList();
};

std::ostream& operator<<(std::ostream &out, GraphComponentMap& component_map);

typedef std::shared_ptr<GraphComponentMap> GraphComponentMapPtr;