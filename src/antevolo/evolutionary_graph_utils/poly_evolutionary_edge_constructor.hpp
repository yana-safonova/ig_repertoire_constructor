#pragma once

#include "../antevolo_config.hpp"
#include "base_evolutionary_edge.hpp"
#include "directed_evolutionary_edge.hpp"
#include "undirected_evolutionary_edge.hpp"
#include "intersected_evolutionary_edge.hpp"

namespace antevolo {
    class PolyEvolutionaryEdgeConstructor {
    public:
        virtual std::shared_ptr<BaseEvolutionaryEdge> ConstructEdge(const annotation_utils::AnnotatedClone &src_clone,
                              const annotation_utils::AnnotatedClone &dst_clone,
                                               size_t src_num,
                                               size_t dst_num) const = 0;

        virtual ~PolyEvolutionaryEdgeConstructor() { }
    };

    class PolySimpleEvolutionaryEdgeConstructor : public PolyEvolutionaryEdgeConstructor {
        const AntEvoloConfig::AlgorithmParams::EdgeConstructionParams &params_;

    public:
        PolySimpleEvolutionaryEdgeConstructor(const AntEvoloConfig::AlgorithmParams::EdgeConstructionParams &params) :
                params_(params) { }

        std::shared_ptr<BaseEvolutionaryEdge> ConstructEdge(const annotation_utils::AnnotatedClone &src_clone,
                                       const annotation_utils::AnnotatedClone &dst_clone,
                                       size_t src_num,
                                       size_t dst_num) const;
    };

    class PolyVJEvolutionaryEdgeConstructor : public PolyEvolutionaryEdgeConstructor {
        const AntEvoloConfig::AlgorithmParams::EdgeConstructionParams &params_;

    public:
        PolyVJEvolutionaryEdgeConstructor(const AntEvoloConfig::AlgorithmParams::EdgeConstructionParams &params) :
                params_(params) { }

        std::shared_ptr<BaseEvolutionaryEdge> ConstructEdge(const annotation_utils::AnnotatedClone &src_clone,
                                       const annotation_utils::AnnotatedClone &dst_clone,
                                       size_t src_num,
                                       size_t dst_num) const;
    };
}