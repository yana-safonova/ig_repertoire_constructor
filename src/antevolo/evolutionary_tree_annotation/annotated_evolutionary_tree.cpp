#include <verify.hpp>
#include "annotated_evolutionary_tree.hpp"
#include "tree_based_shm_annotator.hpp"

namespace antevolo {
    void AnnotatedEvolutionaryTree::CheckConsistencyFatal(size_t clone_id) const {
        VERIFY_MSG(tree_ptr_->ContainsClone(clone_id), "Evolutionary tree does not contain clone " << clone_id);
    }

    void AnnotatedEvolutionaryTree::AddSegmentSHMForClone(size_t clone_id, annotation_utils::SHM shm) {
        CheckConsistencyFatal(clone_id);
        if(segment_shms_map_.find(clone_id) == segment_shms_map_.end()) {
            segment_shms_map_[clone_id] = std::vector<EvolutionaryAnnotatedSHM>();
        }
        auto annotated_shm = shm_annotator_.GetAnnotation(clone_id, shm);
        segment_shms_map_[clone_id].push_back(annotated_shm);
        all_segment_shms_.push_back(annotated_shm);
    }

    size_t AnnotatedEvolutionaryTree::NumSynonymousSHMs() const {
        size_t num_synonymous_shms = 0;
        for(auto it = segment_shms_map_.begin(); it != segment_shms_map_.end(); it++) {
            auto shm_vector = it->second;
            for(auto it2 = shm_vector.begin(); it2 != shm_vector.end(); it2++)
                if(it2->synonymous)
                    num_synonymous_shms++;
        }
        for(auto it = cdr3_shms_map_.begin(); it != cdr3_shms_map_.end(); it++) {
            auto shm_vector = it->second;
            for(auto it2 = shm_vector.begin(); it2 != shm_vector.end(); it2++)
                if(it2->IsSynonymous())
                    num_synonymous_shms++;
        }
        return num_synonymous_shms;
    }

    size_t AnnotatedEvolutionaryTree::NumSynonymousWrtGermlineSHMs() const {
        size_t num_synonymous_shms = 0;
        for(auto it = segment_shms_map_.begin(); it != segment_shms_map_.end(); it++) {
            auto shm_vector = it->second;
            for(auto it2 = shm_vector.begin(); it2 != shm_vector.end(); it2++)
                if(it2->shm.IsSynonymous())
                    num_synonymous_shms++;
        }
        return num_synonymous_shms;
    }

    size_t AnnotatedEvolutionaryTree::RootDepth() const {
        size_t root_id = tree_ptr_->GetRoot();
        return (*clone_set_prt_)[root_id].VSHMs().size() + (*clone_set_prt_)[root_id].JSHMs().size();
    }

    size_t AnnotatedEvolutionaryTree::NumAddedSHMs() const {
        return NumUniqueSHms() - RootDepth();
    }

    size_t AnnotatedEvolutionaryTree::SHMDepth() const {
        size_t shm_depth = 0;
        for(auto it = tree_ptr_->c_vertex_begin(); it != tree_ptr_->c_vertex_end(); it++) {
            if(tree_ptr_->IsLeaf(*it)) {
                auto leaf_clone = (*clone_set_prt_)[*it];
                shm_depth = std::max(shm_depth, leaf_clone.VSHMs().size() + leaf_clone.JSHMs().size() +
                        GetNumCDR3SHMsFromCloneToRoot(*it) - RootDepth());
            }
        }
        return shm_depth;
    }

    void AnnotatedEvolutionaryTree::AddCDR3SHMForClone(size_t src_id, size_t dst_id, size_t src_pos, size_t dst_pos) {
        CheckConsistencyFatal(src_id);
        CheckConsistencyFatal(dst_id);
        if(cdr3_shms_map_.find(dst_id) != cdr3_shms_map_.end()) {
            cdr3_shms_map_[dst_id] = std::vector<CDR3SHM>();
        }
        CDR3SHM cdr3_shm(src_pos, dst_pos,
                         (*clone_set_prt_)[src_id].Read().seq[src_pos],
                         (*clone_set_prt_)[dst_id].Read().seq[dst_pos],
                         (*clone_set_prt_)[src_id].GetAminoAcidByNucleotidePos(src_pos),
                         (*clone_set_prt_)[dst_id].GetAminoAcidByNucleotidePos(dst_pos));
        cdr3_shms_map_[dst_id].push_back(cdr3_shm);
        all_cdr3_shms_.push_back(cdr3_shm);
    }

    size_t AnnotatedEvolutionaryTree::GetNumCDR3SHMsForClone(size_t clone_id) const {
        CheckConsistencyFatal(clone_id);
        if(cdr3_shms_map_.find(clone_id) == cdr3_shms_map_.end())
            return 0;
        return cdr3_shms_map_.at(clone_id).size();
    }

    size_t AnnotatedEvolutionaryTree::GetNumCDR3SHMsFromCloneToRoot(size_t clone_id) const {
        CheckConsistencyFatal(clone_id);
        size_t cdr3_depth = 0;
        size_t current_clone = clone_id;
        while(!tree_ptr_->IsRoot(current_clone)) {
            cdr3_depth += GetNumCDR3SHMsForClone(current_clone);
            current_clone = tree_ptr_->GetParentEdge(current_clone).src_clone_num;
        }
        return cdr3_depth;
    }
}