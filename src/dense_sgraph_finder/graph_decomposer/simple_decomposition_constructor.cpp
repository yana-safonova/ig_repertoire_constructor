#include "simple_decomposition_constructor.hpp"

using namespace dense_subgraph_finder;

void SimpleDecompositionConstructor::CreateFirstSet() {
    size_t first_index = permutation_prt_->Reverse()[0];
    decomposition_ptr_->SetClass(first_index, 0);
}

// vertex should be consistent with ids in hamming_graph_ptr
double SimpleDecompositionConstructor::ComputeEdgePercToPreviousSet(size_t vertex) {
    double edges_perc = 0;
    for(size_t i = hamming_graph_ptr_->RowIndex()[vertex];i < hamming_graph_ptr_->RowIndex()[vertex + 1]; i++) {
        size_t old_neigh = hamming_graph_ptr_->Col()[i];
        if(collapsed_struct_->OldIndexValid(old_neigh)) {
            size_t new_neigh = collapsed_struct_->NewVerticesList()[old_neigh];
            if(decomposition_ptr_->LastClassContains(new_neigh))
                edges_perc += 1.0;
        }
    }
    for(size_t i = hamming_graph_ptr_->RowIndexT()[vertex]; i < hamming_graph_ptr_->RowIndexT()[vertex + 1]; i++) {
        size_t old_neigh = hamming_graph_ptr_->ColT()[i];
        if(collapsed_struct_->OldIndexValid(old_neigh)) {
            size_t new_neigh = collapsed_struct_->NewVerticesList()[old_neigh];
            if(decomposition_ptr_->LastClassContains(new_neigh))
                edges_perc += 1.0;
        }
    }
    return edges_perc / double(decomposition_ptr_->LastClassSize());
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
        size_t old_graph_index = collapsed_struct_->OldVerticesList()[old_perm_index];
        TRACE("Old graph index: " << old_graph_index);
        double edge_perc = ComputeEdgePercToPreviousSet(old_graph_index);
        TRACE("Edge %: " << edge_perc);
        assert(edge_perc <= 1.0);
        if(edge_perc >= edge_perc_threshold_) {
            decomposition_ptr_->SetClass(old_perm_index, cur_set);
            TRACE("Set " << cur_set << " was updated");
        }
        else {
            cur_set++;
            decomposition_ptr_->SetClass(old_perm_index, cur_set);
            TRACE("New set was created");
        }
    }
    DEBUG(decomposition_ptr_->Size() << " classes were constructed");
    DEBUG("Maximal class size: " << decomposition_ptr_->MaxClassSize());
    return decomposition_ptr_;
    //ComputeStats();
}