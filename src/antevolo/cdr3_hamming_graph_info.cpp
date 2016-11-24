#include "cdr3_hamming_graph_info.hpp"

namespace antevolo {

    size_t CDR3HammingGraphInfo::GetOldIndexByNewIndex(size_t i) {
        return graph_component_map_.GetOldVertexByNewVertex(component_id_, i);
    }
    size_t CDR3HammingGraphInfo::GetNewIndexByOldIndex(size_t old_index) {
        return graph_component_map_.GetNewVertexByOldVertex(old_index);
    }
    const std::vector<size_t>& CDR3HammingGraphInfo::GetClonesByCDR3(std::string cdr3) {
        return cdr3_to_indices_vector_map_.find(cdr3)->second;
    }
    const std::string& CDR3HammingGraphInfo::GetCDR3ByOldIndex(size_t old_index) {
        return unique_cdr3s_[old_index];
    }
    size_t CDR3HammingGraphInfo::GetOldIndexByCDR3(const std::string& cdr3) {
        return cdr3_to_old_index_map_.find(cdr3)->second;
    }
    const std::vector<size_t>& CDR3HammingGraphInfo::GetClonesByOldIndex(size_t old_index) {
        return GetClonesByCDR3(GetCDR3ByOldIndex(old_index));
    }
    SparseGraph::EdgesIterator CDR3HammingGraphInfo::GetSimilarCDR3sBeginByOldIndex(size_t old_index) {
        size_t new_index = graph_component_map_.GetNewVertexByOldVertex(old_index);
        return hg_component_->VertexEdges(new_index).begin();
    }
    SparseGraph::EdgesIterator CDR3HammingGraphInfo::GetSimilarCDR3sEndByOldIndex(size_t old_index) {
        size_t new_index = graph_component_map_.GetNewVertexByOldVertex(old_index);
        return hg_component_->VertexEdges(new_index).end();
    }

}