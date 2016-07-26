#pragma once

#include "annotated_clone.hpp"
#include "aa_annotation/aa_calculator.hpp"
#include "shm_annotation/shm_calculator.hpp"

namespace annotation_utils {
    class AnnotatedCloneCalculator {
        std::shared_ptr<BaseAACalculator> aa_calculator_ptr_;
        std::shared_ptr<BaseSHMCalculator> v_shm_calculator_ptr_;
        std::shared_ptr<BaseSHMCalculator> j_shm_calculator_ptr_;

    public:
        AnnotatedCloneCalculator(std::shared_ptr<BaseAACalculator> aa_calculator_ptr,
                                 std::shared_ptr<BaseSHMCalculator> v_shm_calculator_ptr,
                                 std::shared_ptr<BaseSHMCalculator> j_shm_calculator_ptr) :
                aa_calculator_ptr_(aa_calculator_ptr),
                v_shm_calculator_ptr_(v_shm_calculator_ptr),
                j_shm_calculator_ptr_(j_shm_calculator_ptr) {
        }

        AnnotatedClone ComputeAnnotatedClone(const core::Read &read,
                                             CDRLabeling cdr_labeling,
                                             alignment_utils::ImmuneGeneReadAlignment v_alignment,
                                             alignment_utils::ImmuneGeneReadAlignment j_alignment);
    };
}