#pragma once

#include "base_evolutionary_edge.hpp"

namespace  antevolo {

    class UndirectedEvolutionaryEdge : public BaseEvolutionaryEdge {
        size_t num_added_v_shms;
        size_t num_intersected_v_shms;
        size_t num_added_j_shms;
        size_t num_intersected_j_shms;
        size_t num_added_shms;
        size_t num_intersected_shms;

    public:
        UndirectedEvolutionaryEdge(const annotation_utils::AnnotatedClone &src_clone,
                                   const annotation_utils::AnnotatedClone &dst_clone,
                                   size_t src_num, size_t dst_num) : BaseEvolutionaryEdge(src_clone, dst_clone,
                                                                          src_num, dst_num) {
            edge_type = EvolutionaryEdgeType::UndirectedEdgeType;
            cdr3_distance = HammingDistance(this->src_clone->CDR3(), this->dst_clone->CDR3());
            num_added_v_shms = 0;
            num_intersected_v_shms = this->src_clone->VSHMs().size();
            num_added_j_shms = 0;
            num_intersected_j_shms = this->src_clone->JSHMs().size();
            num_added_shms = num_added_v_shms + num_added_j_shms;
            num_intersected_shms = num_intersected_v_shms + num_intersected_j_shms;
        }
        size_t Length() const override {
            return cdr3_distance;
        }

        bool IsUndirected() const override { return true; }

        std::string TypeString()  const override { return "undirected"; }

        size_t NumAddedShms() const override { return num_added_shms; }

        size_t NumSharedShms() const override { return num_intersected_shms; }

    };

}