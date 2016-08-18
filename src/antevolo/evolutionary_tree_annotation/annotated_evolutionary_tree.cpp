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
        unique_shms_[clone_id].push_back(shm_annotator_.GetAnnotation(clone_id, shm));
        num_unique_shms_++;
    }
}