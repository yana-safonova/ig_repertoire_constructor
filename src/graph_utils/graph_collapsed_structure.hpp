#pragma once

#include "sparse_graph.hpp"

class GraphCollapsedStructure {
    SparseGraphPtr hamming_graph_ptr_;

    vector<size_t> old_vertices_list_; // old_vertices_list[i] shows old number of new vertex i
    vector<size_t> new_vertices_list_; // new_vertices_list[i] shows new number of old vertex i
    vector<size_t> main_vertices_map_; // main_vertices_map[i] shows index of main vertex (old indices)
    vector<size_t> multiplicity_;

    void InitializeMainVerticesMap();

    void InitializeOldVerticesList();

    void InitializeNewVerticesList();

    void Initialize();

public:
    GraphCollapsedStructure(SparseGraphPtr hamming_graph_ptr) : hamming_graph_ptr_(hamming_graph_ptr) {
        Initialize();
    }

    const vector<size_t>& MainVerticesMap() const { return main_vertices_map_; }

    const vector<size_t>& OldVerticesList() const { return old_vertices_list_; }

    const vector<size_t>& NewVerticesList() const { return new_vertices_list_; }

    size_t NumCollapsedVertices() const { return old_vertices_list_.size(); }

    bool OldIndexValid(size_t index) {
        assert(index < new_vertices_list_.size());
        return new_vertices_list_[index] != size_t(-1);
    }

    bool NewIndexValid(size_t index) {
        return index != size_t(-1);
    }

    size_t GetMultiplicityOf(size_t index) {
        assert(index < multiplicity_.size());
        return multiplicity_[index];
    }

    bool VertexIsMain(size_t index) {
        assert(index < main_vertices_map_.size());
        return index == main_vertices_map_[index];
    }

    size_t GetMainVertexIndex(size_t index) {
        assert(index < main_vertices_map_.size());
        return main_vertices_map_[index];
    }

    size_t NewIndexOfOldVertex(size_t old_vertex) {
        assert(old_vertex < main_vertices_map_.size());
        return new_vertices_list_[main_vertices_map_[old_vertex]];
    }

    size_t OldIndexOfNewVertex(size_t new_vertex) {
        assert(new_vertex < old_vertices_list_.size());
        return old_vertices_list_[new_vertex];
    }

    size_t NumberNewVertices() { return old_vertices_list_.size(); }

    size_t NumberOldVertices() { return main_vertices_map_.size(); }

    size_t MultiplicityOfNewVertex(size_t new_index) {
        assert(new_index < old_vertices_list_.size());
        return multiplicity_[old_vertices_list_[new_index]];
    }

    size_t MultiplicityOfOldVertex(size_t old_vertex) {
        assert(old_vertex < multiplicity_.size());
        return multiplicity_[old_vertex];
    }

    size_t NumberCollapsedEdges(SparseGraphPtr hamming_graph_ptr);
};

typedef shared_ptr<GraphCollapsedStructure> GraphCollapsedStructurePtr;

ostream& operator<<(ostream& out, GraphCollapsedStructure collapsed_structure);
