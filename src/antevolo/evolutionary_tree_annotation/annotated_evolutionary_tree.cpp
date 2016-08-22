#include <verify.hpp>
#include "annotated_evolutionary_tree.hpp"
#include "tree_based_shm_annotator.hpp"

namespace antevolo {
    void AnnotatedEvolutionaryTree::CheckConsistencyFatal(size_t clone_id) {
        VERIFY_MSG(tree_ptr_->ContainsClone(clone_id), "Evolutionary tree does not contain clone " << clone_id);
    }

    void AnnotatedEvolutionaryTree::AddSHMForClone(size_t clone_id, annotation_utils::SHM shm) {
        CheckConsistencyFatal(clone_id);
        if(unique_shms_.find(clone_id) == unique_shms_.end()) {
            unique_shms_[clone_id] = std::vector<EvolutionaryAnnotatedSHM>();
        }
        auto annotated_shm = shm_annotator_.GetAnnotation(clone_id, shm);
        unique_shms_[clone_id].push_back(annotated_shm);
        all_unique_shms_.push_back(annotated_shm);
    }

    size_t AnnotatedEvolutionaryTree::NumSynonymousSHMs() const {
        size_t num_synonymous_shms = 0;
        for(auto it = unique_shms_.begin(); it != unique_shms_.end(); it++) {
            auto shm_vector = it->second;
            for(auto it2 = shm_vector.begin(); it2 != shm_vector.end(); it2++)
                if(it2->synonymous)
                    num_synonymous_shms++;
        }
        return num_synonymous_shms;
    }

    size_t AnnotatedEvolutionaryTree::NumSynonymousWrtGermlineSHMs() const {
        size_t num_synonymous_shms = 0;
        for(auto it = unique_shms_.begin(); it != unique_shms_.end(); it++) {
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
}