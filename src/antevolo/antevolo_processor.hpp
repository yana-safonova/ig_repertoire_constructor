#pragma once

#include "evolutionary_tree_storage.hpp"
#include "antevolo_config.hpp"
#include <annotation_utils/annotated_clone_set.hpp>

namespace antevolo {
    class AntEvoloProcessor {
        const AntEvoloConfig& config_;
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;

        std::vector<EvolutionaryTreeStorage> thread_tree_storages_;

        EvolutionaryTreeStorage JoinEvolutionaryStoragesFromThreads();

    public:
        AntEvoloProcessor(const AntEvoloConfig& config,
                          const annotation_utils::CDRAnnotatedCloneSet &clone_set) :
                config_(config),
                clone_set_(clone_set) {
            for(int i = 0; i < config.run_params.num_threads; i++)
                thread_tree_storages_.push_back(EvolutionaryTreeStorage(clone_set));
        }

        EvolutionaryTreeStorage ConstructClonalTrees();
    };
}