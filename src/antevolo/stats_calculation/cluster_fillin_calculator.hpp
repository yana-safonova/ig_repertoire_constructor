#pragma once

#include <annotation_utils/annotated_clone_set.hpp>
#include <boost/unordered/unordered_set.hpp>
#include <evolutionary_graph_utils/evolutionary_edge_constructor.hpp>
#include "antevolo_config.hpp"

namespace antevolo {

    class ClusterFillinCalculator {
        const annotation_utils::CDRAnnotatedCloneSet& clone_set_;
        const AntEvoloConfig& config_;

    public:
        ClusterFillinCalculator(const annotation_utils::CDRAnnotatedCloneSet &clone_set,
                                const AntEvoloConfig& config) :
                clone_set_(clone_set),
                config_(config){}

        size_t CalculateE(boost::unordered_set<size_t>);

    private:
        std::shared_ptr<EvolutionaryEdgeConstructor> GetEdgeConstructor() {
            EvolutionaryEdgeConstructor* ptr = new VJEvolutionaryEdgeConstructor(
                    config_.algorithm_params.edge_construction_params);
            return std::shared_ptr<EvolutionaryEdgeConstructor>(ptr);
        }
    };
}

