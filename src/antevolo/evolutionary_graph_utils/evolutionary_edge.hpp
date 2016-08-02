#pragma once

#include <annotation_utils/annotated_clone.hpp>

namespace antevolo {
    enum EvolutionaryEdgeType { UnknownEdgeType, DirectedEdgeType, UndirectedEdgeType, IntersectedEdgeType };

    std::ostream& operator<<(std::ostream &out, const EvolutionaryEdgeType& edge_type);

    struct EvolutionaryEdge {
        EvolutionaryEdgeType edge_type;
        std::shared_ptr<annotation_utils::AnnotatedClone> src_clone;
        std::shared_ptr<annotation_utils::AnnotatedClone> dst_clone;
        size_t src_clone_num;
        size_t dst_clone_num;

        size_t num_added_v_shms;
        size_t num_intersected_v_shms;
        size_t num_added_j_shms;
        size_t num_intersected_j_shms;
        size_t num_added_shms;
        size_t num_intersected_shms;
        size_t cdr3_distance;
        size_t weight;

        void InitializeFields();

    public:
        EvolutionaryEdge(EvolutionaryEdgeType edge_type,
                         const annotation_utils::AnnotatedClone &src_clone,
                         const annotation_utils::AnnotatedClone &dst_clone,
                         size_t src_num,
                         size_t dst_num,
                         size_t intersected_edge_coeff) : edge_type(edge_type),
                                                          src_clone(std::make_shared<annotation_utils::AnnotatedClone>(src_clone)),
                                                          dst_clone(std::make_shared<annotation_utils::AnnotatedClone>(dst_clone)),
                                                          src_clone_num(src_num),
                                                          dst_clone_num(dst_num),
                                                          intersected_edge_coeff_(intersected_edge_coeff) {
            InitializeFields();
        }

        EvolutionaryEdge() {}

        EvolutionaryEdge(const EvolutionaryEdge& edge) {
            edge_type = edge.edge_type;
            src_clone = edge.src_clone;
            dst_clone = edge.dst_clone;
            src_clone_num = edge.src_clone_num;
            dst_clone_num = edge.dst_clone_num;
            num_added_v_shms = edge.num_added_v_shms;
            num_intersected_v_shms = edge.num_intersected_v_shms;
            num_added_j_shms = edge.num_added_j_shms;
            num_intersected_j_shms = edge.num_intersected_j_shms;
            num_added_shms = edge.num_added_shms;
            num_intersected_shms = edge.num_intersected_shms;
            cdr3_distance = edge.cdr3_distance;
            weight = edge.weight;
        }

        bool Empty() const { return edge_type == EvolutionaryEdgeType::UnknownEdgeType; }

        bool IsDirected() const { return edge_type == EvolutionaryEdgeType::DirectedEdgeType; };

        bool IsUndirected() const { return edge_type == EvolutionaryEdgeType::UndirectedEdgeType; };

        bool IsSynonymous() const;

        bool operator==(const EvolutionaryEdge &edge) const {
            return src_clone == edge.src_clone and dst_clone == edge.dst_clone and edge_type == edge.edge_type;
        }

    private:
        size_t intersected_edge_coeff_;
    };

    class WeightEvolutionaryEdgeComparator {
    public:
        bool operator() (const EvolutionaryEdge& e1, const EvolutionaryEdge& e2) {
            return e1.weight < e2.weight;
        }
    };

    std::ostream& operator<<(std::ostream& out, const EvolutionaryEdge &edge);
}