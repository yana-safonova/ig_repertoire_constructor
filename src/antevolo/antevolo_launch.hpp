#pragma once

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <shm_model_utils/shm_model_edge_weight_calculator.hpp>
#include "antevolo_config.hpp"
#include "annotation_utils/annotated_clone_set.hpp"
#include "annotated_clone_by_read_constructor.hpp"
#include "evolutionary_tree_storage.hpp"

namespace antevolo {
    class AntEvoloLaunch {
        const AntEvoloConfig& config_;

        ShmModelEdgeWeightCalculator ShmModelPosteriorCalculation(
                const annotation_utils::AnnotatedCloneSet<annotation_utils::AnnotatedClone>&);
        std::vector<boost::unordered_set<size_t>> ReadClusters(
                const boost::unordered_map<std::string, size_t>& read_name_to_index);

        void LaunchDefault(const AnnotatedCloneByReadConstructor& clone_by_read_constructor,
                           const annotation_utils::CDRAnnotatedCloneSet& annotated_clone_set,
                           size_t total_number_of_reads);

        void LaunchEvoQuast(const annotation_utils::CDRAnnotatedCloneSet& clone_set);

        void AnalyzeParallelEvolution(const EvolutionaryTreeStorage& trees);

    public:
        AntEvoloLaunch(const AntEvoloConfig& config) : config_(config) { }

        void Launch();
    };
}