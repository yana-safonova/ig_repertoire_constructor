#pragma once

#include "evolutionary_tree_storage.hpp"
#include "antevolo_config.hpp"
#include <annotation_utils/annotated_clone_set.hpp>

namespace antevolo {
    class AntEvoloProcessor {
        const AntEvoloConfig& config_;
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;

        std::vector<EvolutionaryTreeStorage> thread_tree_storages_;

        // todo: refactor this methods: too many arguments
        std::string GetTreeOutputFname(std::string output_dir, size_t index1, size_t index2, size_t v_num, size_t e_num);
        std::string GetTreeClonesOutputFname(std::string output_dir, size_t index1, size_t index2, size_t v_num, size_t e_num);

        EvolutionaryTreeStorage JoinEvolutionaryStoragesFromThreads();

    public:
        AntEvoloProcessor(const AntEvoloConfig& config,
                          const annotation_utils::CDRAnnotatedCloneSet &clone_set) :
                config_(config),
                clone_set_(clone_set) {
            for(size_t i = 0; i < config.run_params.num_threads; i++)
                thread_tree_storages_.push_back(EvolutionaryTreeStorage(clone_set));
        }

        EvolutionaryTreeStorage ConstructClonalTrees();
    };
}