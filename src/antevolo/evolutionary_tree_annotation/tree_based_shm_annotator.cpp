#include "tree_based_shm_annotator.hpp"

namespace antevolo {
    EvolutionaryAnnotatedSHM TreeBasedSHMAnnotator::GetAnnotationForRoot(size_t root_id, annotation_utils::SHM shm) {
        VERIFY_MSG(tree_ptr_->IsRoot(root_id), "Vertex " << root_id << " is not root");
        return EvolutionaryAnnotatedSHM(shm, false, tree_ptr_->NumVertices() > 1);
    }

    // todo: write unit test
    bool TreeBasedSHMAnnotator::SHMIsSynonymousWrtToParent(size_t dst_clone, annotation_utils::SHM shm) {
        size_t src_clone = tree_ptr_->GetParentEdge(dst_clone)->SrcNum();
        auto src_gene_alignment = (*clone_set_ptr_)[src_clone].GetAlignmentBySegment(shm.segment_type);
        size_t src_nucl_pos = src_gene_alignment.QueryPositionBySubjectPosition(shm.gene_nucl_pos);
        return shm.read_aa == (*clone_set_ptr_)[src_clone].GetAminoAcidByNucleotidePos(src_nucl_pos);
    }

    EvolutionaryAnnotatedSHM TreeBasedSHMAnnotator::GetAnnotationForNodeWithParent(size_t clone_id,
                                                                                   annotation_utils::SHM shm) {
        VERIFY_MSG(tree_ptr_->ContainsClone(clone_id), "Tree does not contain clone " << clone_id);
        VERIFY_MSG(!tree_ptr_->IsRoot(clone_id), "Node " << clone_id << " is root of the tree");
        return EvolutionaryAnnotatedSHM(shm, SHMIsSynonymousWrtToParent(clone_id, shm), false); //tree_ptr_->IsLeaf(clone_id));
    }

    EvolutionaryAnnotatedSHM TreeBasedSHMAnnotator::GetAnnotation(size_t clone_id, annotation_utils::SHM shm) {
        VERIFY_MSG(tree_ptr_->ContainsClone(clone_id), "Tree does not contain clone " << clone_id);
        if(tree_ptr_->IsRoot(clone_id))
            return GetAnnotationForRoot(clone_id, shm);
        return GetAnnotationForNodeWithParent(clone_id, shm);
    }
}