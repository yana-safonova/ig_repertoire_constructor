#pragma once

#include "../antevolo_config.hpp"
#include "evolutionary_edge.hpp"

namespace antevolo {
    class EvolutionaryEdgeConstructor {
    public:
        virtual EvolutionaryEdge ConstructEdge(const annotation_utils::AnnotatedClone &src_clone,
                              const annotation_utils::AnnotatedClone &dst_clone,
                                               size_t src_num,
                                               size_t dst_num) const = 0;

        virtual ~EvolutionaryEdgeConstructor() { }
    };

    class SimpleEvolutionaryEdgeConstructor : public EvolutionaryEdgeConstructor {
        const AntEvoloConfig::AlgorithmParams::EdgeConstructionParams &params_;

    public:
        SimpleEvolutionaryEdgeConstructor(const AntEvoloConfig::AlgorithmParams::EdgeConstructionParams &params) :
                params_(params) { }

        EvolutionaryEdge ConstructEdge(const annotation_utils::AnnotatedClone &src_clone,
                                       const annotation_utils::AnnotatedClone &dst_clone,
                                       size_t src_num,
                                       size_t dst_num) const;
    };
}