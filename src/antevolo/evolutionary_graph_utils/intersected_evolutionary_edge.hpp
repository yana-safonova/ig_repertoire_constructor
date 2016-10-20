#pragma once

#include "base_evolutionary_edge.hpp"

namespace  antevolo {

    class IntersectedEvolutionaryEdge : public BaseEvolutionaryEdge {
        size_t num_individual_v_shms;
        size_t num_intersected_v_shms;
        size_t num_individual_j_shms;
        size_t num_intersected_j_shms;
        size_t num_individual_shms;
        size_t num_intersected_shms;

    public:
        IntersectedEvolutionaryEdge(const annotation_utils::AnnotatedClone &src_clone,
                                    const annotation_utils::AnnotatedClone &dst_clone,
                                    size_t src_num, size_t dst_num,
                                    size_t intersected_edge_coeff) : BaseEvolutionaryEdge(edge_type,
                                                                          src_clone,
                                                                          dst_clone,
                                                                          src_num,
                                                                          dst_num) {
            edge_type = EvolutionaryEdgeType::IntersectedEdgeType;
            cdr3_distance = HammingDistance(this->src_clone->CDR3(), this->dst_clone->CDR3());
            //V
            num_intersected_v_shms = annotation_utils::SHMComparator::GetNumberOfIntersections(this-src_clone->VSHMs(),
                                                                                               this-dst_clone->VSHMs());
            num_individual_v_shms = this-src_clone->VSHMs().size() + this-dst_clone->VSHMs().size() - 2 * num_intersected_v_shms;

            //J
            num_intersected_j_shms = annotation_utils::SHMComparator::GetNumberOfIntersections(this-src_clone->JSHMs(),
                                                                                               this-dst_clone->JSHMs());
            num_individual_j_shms = this-src_clone->JSHMs().size() + this-dst_clone->JSHMs().size() - 2 * num_intersected_j_shms;

            //sum
            num_individual_shms = num_individual_v_shms + num_individual_v_shms;
            num_intersected_shms = num_intersected_v_shms + num_intersected_j_shms;
        }
        size_t Length() const {
            return cdr3_distance+num_individual_shms;
        }

        string TypeString() const { return "intersected"; }

        size_t NumAddedShms() const { return num_individual_shms; }

        size_t NumSharedShms() const { return num_intersected_shms; }

    };

}