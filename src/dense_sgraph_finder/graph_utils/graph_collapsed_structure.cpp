#include "graph_collapsed_structure.hpp"

void GraphCollapsedStructure::InitializeMainVerticesMap() {
    vector<size_t> main_vertices_tree;
    for(size_t i = 0; i < hamming_graph_ptr_->N(); i++) {
        multiplicity_.push_back(1);
        main_vertices_tree.push_back(i);
    }

    for(size_t i = 0; i < hamming_graph_ptr_->N(); i++)
        for(size_t j = hamming_graph_ptr_->RowIndex()[i]; j < hamming_graph_ptr_->RowIndex()[i + 1];j++)
            if(hamming_graph_ptr_->Dist()[j] == 0)
                main_vertices_tree[hamming_graph_ptr_->Col()[j]] = i;

    for(size_t i = 0; i < main_vertices_tree.size(); i++) {
        size_t cur_vertex = i;
        while(cur_vertex != main_vertices_tree[cur_vertex])
            cur_vertex = main_vertices_tree[cur_vertex];
        main_vertices_map_.push_back(cur_vertex);
        multiplicity_[cur_vertex]++;
    }
}

void GraphCollapsedStructure::InitializeOldVerticesList() {
    for(size_t i = 0; i < main_vertices_map_.size(); i++)
        if(main_vertices_map_[i] == i)
            old_vertices_list_.push_back(i);
}

void GraphCollapsedStructure::InitializeNewVerticesList() {
    for(size_t i = 0; i < hamming_graph_ptr_->N(); i++)
        new_vertices_list_.push_back(size_t(-1));
    for(size_t i = 0; i < old_vertices_list_.size(); i++)
        new_vertices_list_[old_vertices_list_[i]] = i;
}

void GraphCollapsedStructure::Initialize() {
    InitializeMainVerticesMap();
    InitializeOldVerticesList();
    InitializeNewVerticesList();
}

size_t GraphCollapsedStructure::NumberCollapsedEdges(SparseGraphPtr hamming_graph_ptr) {
    size_t collapsed_edges = 0;
    for(size_t i = 0; i < hamming_graph_ptr->N(); i++)
        for(size_t j = hamming_graph_ptr->RowIndex()[i]; j < hamming_graph_ptr->RowIndex()[i + 1]; j++) {
            size_t v1 = i;
            size_t v2 = hamming_graph_ptr->Col()[j];
            if(VertexIsMain(v1) and VertexIsMain(v2)) {
                collapsed_edges++;
            }
        }
    return collapsed_edges;
}

ostream& operator<<(ostream& out, GraphCollapsedStructure collapsed_structure) {
    out << "Old vertex - new vertex" << endl;
    for(size_t i = 0; i < collapsed_structure.NumberNewVertices(); i++)
        out << collapsed_structure.OldVerticesList()[i] << "\t" << i << endl;
    out << "Old vertex - main vertex" << endl;
    for(size_t i = 0; i < collapsed_structure.NumberOldVertices(); i++)
        out << i << "\t" << collapsed_structure.GetMainVertexIndex(i) << endl;
    return out;
}