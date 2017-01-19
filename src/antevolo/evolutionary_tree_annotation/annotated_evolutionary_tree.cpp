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
        unique_segment_shms_.insert(annotated_shm);
    }

    size_t AnnotatedEvolutionaryTree::NumSynonymousSegmentSHMs() const {
        size_t num_synonymous_shms = 0;
        for(auto it = unique_segment_shms_.begin(); it != unique_segment_shms_.end(); it++) {
            if(it->synonymous)
                num_synonymous_shms++;
        }
        return num_synonymous_shms;
    }

    size_t AnnotatedEvolutionaryTree::NumSynonymousCDR3SHMs() const {
        size_t num_synonymous_shms = 0;
        for(auto it = unique_cdr3_shms_.begin(); it != unique_cdr3_shms_.end(); it++) {
            if (it->IsSynonymous())
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
        const auto& clone_set = GetCloneSet();
        size_t root_id = tree_ptr_->GetRoot();
        return clone_set[root_id].VSHMs().size() + clone_set[root_id].JSHMs().size();
    }

    size_t AnnotatedEvolutionaryTree::NumAddedSHMs() const {
        return NumUniqueSHms() - RootDepth();
    }

    size_t AnnotatedEvolutionaryTree::SHMDepth() const {
        const auto& clone_set = GetCloneSet();
        size_t shm_depth = 0;
        for(auto it = tree_ptr_->c_vertex_begin(); it != tree_ptr_->c_vertex_end(); it++) {
            if(tree_ptr_->IsLeaf(*it)) {
                auto leaf_clone = clone_set[*it];
                shm_depth = std::max(shm_depth, leaf_clone.VSHMs().size() + leaf_clone.JSHMs().size() +
                        GetNumCDR3SHMsFromCloneToRoot(*it) - RootDepth());
            }
        }
        return shm_depth;
    }

    void AnnotatedEvolutionaryTree::AddCDR3SHMForClone(size_t src_id, size_t dst_id, size_t rel_src_pos, size_t rel_dst_pos) {
        const auto& clone_set = GetCloneSet();
        CheckConsistencyFatal(src_id);
        CheckConsistencyFatal(dst_id);
        if(cdr3_shms_map_.find(dst_id) != cdr3_shms_map_.end()) {
            cdr3_shms_map_[dst_id] = std::vector<CDR3SHM>();
        }
        size_t src_pos = clone_set[src_id].CDR3Range().start_pos + rel_src_pos;
        size_t dst_pos = clone_set[dst_id].CDR3Range().start_pos + rel_dst_pos;
        CDR3SHM cdr3_shm(src_pos, dst_pos,
                         rel_src_pos, rel_dst_pos,
                         clone_set[src_id].Read().seq[src_pos],
                         clone_set[dst_id].Read().seq[dst_pos],
                         clone_set[src_id].GetAminoAcidByNucleotidePos(src_pos),
                         clone_set[dst_id].GetAminoAcidByNucleotidePos(dst_pos));
        cdr3_shms_map_[dst_id].push_back(cdr3_shm);
        all_cdr3_shms_.push_back(cdr3_shm);
        unique_cdr3_shms_.insert(cdr3_shm);
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
            current_clone = tree_ptr_->GetParentEdge(current_clone)->SrcNum();
        }
        return cdr3_depth;
    }
}