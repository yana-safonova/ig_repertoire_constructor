#include "annotated_clone_calculator.hpp"

namespace annotation_utils {
    AnnotatedClone AnnotatedCloneCalculator::ComputeAnnotatedClone(const core::Read &read, CDRLabeling cdr_labeling,
                                                                   alignment_utils::ImmuneGeneReadAlignment v_alignment,
                                                                   alignment_utils::ImmuneGeneReadAlignment j_alignment) {
        auto aa_annotation = aa_calculator_ptr_->ComputeAminoAcidAnnotation(read, cdr_labeling);
        AnnotatedClone res(read,
                              cdr_labeling,
                              v_alignment,
                              j_alignment,
                              aa_annotation,
                              v_shm_calculator_ptr_->ComputeSHMs(v_alignment, aa_annotation, cdr_labeling),
                              j_shm_calculator_ptr_->ComputeSHMs(j_alignment, aa_annotation, cdr_labeling));
        return res;
    }
}