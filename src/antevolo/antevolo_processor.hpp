#pragma once

#include "evolutionary_tree_storage.hpp"
#include "antevolo_config.hpp"
#include "clone_set_with_fakes.hpp"
namespace antevolo {
    class AntEvoloProcessor {
        const AntEvoloConfig& config_;
        CloneSetWithFakes& &clone_set_;

        std::vector<EvolutionaryTreeStorage> thread_tree_storages_;

        EvolutionaryTreeStorage JoinEvolutionaryStoragesFromThreads();

    public:
        AntEvoloProcessor(const AntEvoloConfig& config,
                          CloneSetWithFakes& clone_set) :
                config_(config),
                clone_set_(clone_set) {
            for(int i = 0; i < config.run_params.num_threads; i++)
                thread_tree_storages_.push_back(EvolutionaryTreeStorage(clone_set));
        }

        EvolutionaryTreeStorage ConstructClonalTrees();
    };
}