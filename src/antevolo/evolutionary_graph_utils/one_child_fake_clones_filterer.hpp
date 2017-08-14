#pragma once

#include "evolutionary_tree.hpp"

namespace  antevolo {
    class OneChildFakeClonesFilterer {
        const AntEvoloConfig::AlgorithmParams::EdgeConstructionParams params_;
    public:
        OneChildFakeClonesFilterer(const AntEvoloConfig::AlgorithmParams::EdgeConstructionParams params) :
                params_(params) {}
        EvolutionaryTree FilterOneChildFakes(const EvolutionaryTree& tree) const;

        std::shared_ptr<EvolutionaryEdgeConstructor> GetEdgeConstructor() const {
            EvolutionaryEdgeConstructor* ptr = new VJEvolutionaryEdgeConstructor(params_);
            return std::shared_ptr<EvolutionaryEdgeConstructor>(ptr);
        }
    };

}