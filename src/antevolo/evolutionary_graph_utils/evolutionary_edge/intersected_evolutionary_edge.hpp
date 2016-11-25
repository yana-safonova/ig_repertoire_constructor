#pragma once

#include "base_evolutionary_edge.hpp"
#include <annotation_utils/shm_comparator.hpp>

namespace  antevolo {

    class IntersectedEvolutionaryEdge : public BaseEvolutionaryEdge {
        size_t num_individual_v_shms_;
        size_t num_intersected_v_shms_;
        size_t num_individual_j_shms_;
        size_t num_intersected_j_shms_;
        size_t num_individual_shms_;
        size_t num_intersected_shms_;

    public:
        IntersectedEvolutionaryEdge(const annotation_utils::AnnotatedClone &src_clone,
                                    const annotation_utils::AnnotatedClone &dst_clone,
                                    size_t src_num, size_t dst_num,
                                    size_t num_individual_v_shms, size_t num_intersected_v_shms,
                                    size_t num_individual_j_shms, size_t num_intersected_j_shms)
                : BaseEvolutionaryEdge(src_clone,
                                       dst_clone,
                                       src_num,
                                       dst_num),
                  num_individual_v_shms_(num_individual_v_shms),
                  num_intersected_v_shms_(num_intersected_v_shms),
                  num_individual_j_shms_(num_individual_j_shms),
                  num_intersected_j_shms_(num_intersected_j_shms) {
            edge_type = EvolutionaryEdgeType ::IntersectedEdgeType;
            //sum
            num_individual_shms_ = num_individual_v_shms_ + num_individual_v_shms_;
            num_intersected_shms_ = num_intersected_v_shms_ + num_intersected_j_shms_;
        }
        size_t Length() const override {
            return cdr3_distance+num_individual_shms_;
        }

        bool IsIntersected() const override { return true; }

        std::string TypeString() const override { return "intersected"; }

        size_t NumAddedShms() const override { return num_individual_shms_; }

        size_t NumSharedShms() const override { return num_intersected_shms_; }

    };

}