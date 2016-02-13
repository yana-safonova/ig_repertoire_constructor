#include "simple_decomposition_constructor.hpp"

using namespace dense_subgraph_finder;

bool SimpleDecompositionConstructor::VertexIsSupernode(size_t vertex) {
    return graph_ptr_->WeightOfVertex(vertex) >= min_supernode_size_;
}

void SimpleDecompositionConstructor::CreateFirstSet() {
    size_t first_index = permutation_prt_->Reverse()[0];
    decomposition_ptr_->SetClass(first_index, 0);
    current_set_has_snode_ = VertexIsSupernode(first_index);
}

// vertex should be consistent with ids in hamming_graph_ptr
double SimpleDecompositionConstructor::ComputeEdgePercToPreviousSet(size_t vertex) {
    double edges_perc = 0;
    for(size_t i = graph_ptr_->RowIndex()[vertex];i < graph_ptr_->RowIndex()[vertex + 1]; i++) {
        size_t old_neigh = graph_ptr_->Col()[i];
            if(decomposition_ptr_->LastClassContains(old_neigh))
                edges_perc += 1.0;
    }
    for(size_t i = graph_ptr_->RowIndexT()[vertex]; i < graph_ptr_->RowIndexT()[vertex + 1]; i++) {
        size_t old_neigh = graph_ptr_->ColT()[i];
            if(decomposition_ptr_->LastClassContains(old_neigh))
                edges_perc += 1.0;
    }
    return edges_perc / double(decomposition_ptr_->LastClassSize());
}

bool SimpleDecompositionConstructor::GlueVertexWithPreviousSet(size_t vertex) {
    // if vertex is supernode and class already contains supernode, do not glue it
    if(graph_ptr_->WeightOfVertex(vertex) >= min_supernode_size_ and current_set_has_snode_)
        return false;
    double edge_perc = ComputeEdgePercToPreviousSet(vertex);
    TRACE("Edge %: " << edge_perc);
    assert(edge_perc <= 1.0);
    return edge_perc >= edge_perc_threshold_;
}

DecompositionPtr SimpleDecompositionConstructor::CreateDecomposition() {
    DEBUG("Edge % threshold: " << edge_perc_threshold_);
    CreateFirstSet();
    TRACE(*permutation_prt_);
    size_t cur_set = 0;
    for(size_t i = 1; i < permutation_prt_->Size(); i++) {
        TRACE("Index: " << i);
        size_t old_perm_index = permutation_prt_->Reverse()[i];
        TRACE("Old perm index: " << old_perm_index);
        if(GlueVertexWithPreviousSet(old_perm_index)) {
            decomposition_ptr_->SetClass(old_perm_index, cur_set);
            TRACE("Set " << cur_set << " was updated");
            // if newly appended vertex is a supernode, set flag in true
            if(VertexIsSupernode(old_perm_index))
                current_set_has_snode_ = true;
        }
        else {
            cur_set++;
            decomposition_ptr_->SetClass(old_perm_index, cur_set);
            TRACE("New set was created");
            // current class contains the only node: old_perm_index
            // if it is a supernode, set flag in true
            // otherwise in false
            current_set_has_snode_ = VertexIsSupernode(old_perm_index);
        }
    }
    DEBUG(decomposition_ptr_->Size() << " classes were constructed");
    DEBUG("Maximal class size: " << decomposition_ptr_->MaxClassSize());
    return decomposition_ptr_;
}