#pragma once

#include "../antevolo_config.hpp"
#include "evolutionary_graph_utils/evolutionary_edge/base_evolutionary_edge.hpp"
#include "evolutionary_graph_utils/evolutionary_edge/directed_evolutionary_edge.hpp"
#include "evolutionary_graph_utils/evolutionary_edge/undirected_evolutionary_edge.hpp"
#include "evolutionary_graph_utils/evolutionary_edge/intersected_evolutionary_edge.hpp"
#include "evolutionary_graph_utils/evolutionary_edge/reverse_directed_evolutionary_edge.hpp"

namespace antevolo {
    class EvolutionaryEdgeConstructor {
    public:
        virtual std::shared_ptr<BaseEvolutionaryEdge> ConstructEdge(const annotation_utils::AnnotatedClone &src_clone,
                                                                    const annotation_utils::AnnotatedClone &dst_clone,
                                                                    size_t src_num,
                                                                    size_t dst_num) const = 0;

        virtual ~EvolutionaryEdgeConstructor() { }
    };

//    class SimpleEvolutionaryEdgeConstructor : public EvolutionaryEdgeConstructor {
//        const AntEvoloConfig::AlgorithmParams::EdgeConstructionParams &params_;
//
//    public:
//        SimpleEvolutionaryEdgeConstructor(const AntEvoloConfig::AlgorithmParams::EdgeConstructionParams &params) :
//                params_(params) { }
//
//        std::shared_ptr<BaseEvolutionaryEdge> ConstructEdge(const annotation_utils::AnnotatedClone &src_clone,
//                                       const annotation_utils::AnnotatedClone &dst_clone,
//                                       size_t src_num,
//                                       size_t dst_num) const;
//    };

    class VJEvolutionaryEdgeConstructor : public EvolutionaryEdgeConstructor {
        const AntEvoloConfig::AlgorithmParams::EdgeConstructionParams &params_;

    public:
        VJEvolutionaryEdgeConstructor(const AntEvoloConfig::AlgorithmParams::EdgeConstructionParams &params) :
                params_(params) { }

        std::shared_ptr<BaseEvolutionaryEdge> ConstructEdge(const annotation_utils::AnnotatedClone &src_clone,
                                       const annotation_utils::AnnotatedClone &dst_clone,
                                       size_t src_num,
                                       size_t dst_num) const;
    };
}