#include "cdr3_hamming_graph_info.hpp"
#include <log.hpp>
#include <verify.hpp>

namespace antevolo {

    size_t CDR3HammingGraphInfo::GetOldIndexByNewIndex(size_t i) {
        return graph_component_map_.GetOldVertexByNewVertex(component_id_, i);
    }
    size_t CDR3HammingGraphInfo::GetNewIndexByOldIndex(size_t old_index) {
        return graph_component_map_.GetNewVertexByOldVertex(old_index);
    }
    const std::vector<size_t>& CDR3HammingGraphInfo::GetClonesByCDR3(std::string cdr3) {
        //cdr3 = cdr3.substr(0, cdr3.size()-8);
        return cdr3_to_indices_vector_map_.find(cdr3)->second;
    }
    const std::string& CDR3HammingGraphInfo::GetCDR3ByOldIndex(size_t old_index) {
        return unique_cdr3s_[old_index];
    }
    size_t CDR3HammingGraphInfo::GetOldIndexByCDR3(const std::string& cdr3) {
        VERIFY_MSG(cdr3_to_old_index_map_.find(cdr3) != cdr3_to_old_index_map_.end(),
                   "CDR3HammingGraphInfo: failed to find cdr3");
        return cdr3_to_old_index_map_.find(cdr3)->second;
    }

    size_t CDR3HammingGraphInfo::GetOldIndexByCDR3JNucl(std::string cdr3JNucl) {
        VERIFY_MSG(cdr3_to_old_index_map_.find(cdr3JNucl) != cdr3_to_old_index_map_.end(),
                   "CDR3HammingGraphInfo: failed to find cdr3JNucl");
        return cdr3_to_old_index_map_.find(cdr3JNucl)->second;
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
    boost::unordered_set<size_t> CDR3HammingGraphInfo::GetAllClones() {
        boost::unordered_set<size_t> vertices_nums;
        for (size_t i = 0; i < hg_component_->N(); i++) {
//            size_t old_index = graph_component_map_.GetOldVertexByNewVertex(component_id, i);
            size_t old_index = GetOldIndexByNewIndex(i);
            //auto clones_sharing_cdr3 = unique_cdr3s_map_[unique_cdr3Jnucleotides_[old_index]];
            auto clones_sharing_cdr3 = GetClonesByOldIndex(old_index);
            for (size_t clone_num : clones_sharing_cdr3) {
                vertices_nums.insert(clone_num);
            }
        }
        return vertices_nums;
    }



}