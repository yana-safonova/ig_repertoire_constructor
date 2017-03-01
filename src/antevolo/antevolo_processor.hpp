#pragma once

#include "evolutionary_tree_storage.hpp"
#include "antevolo_config.hpp"
#include "clone_set_with_fakes.hpp"
#include "annotated_clone_by_read_constructor.hpp"

namespace antevolo {
    class AntEvoloProcessor {
        const AntEvoloConfig& config_;
        const annotation_utils::CDRAnnotatedCloneSet& clone_set_;
        const AnnotatedCloneByReadConstructor& clone_by_read_constructor_;
        const size_t total_number_of_reads_;

        std::vector<EvolutionaryTreeStorage> thread_tree_storages_;

        EvolutionaryTreeStorage JoinEvolutionaryStoragesFromThreads();

    public:
        AntEvoloProcessor(const AntEvoloConfig& config,
                          const annotation_utils::CDRAnnotatedCloneSet& clone_set,
                          const AnnotatedCloneByReadConstructor& clone_by_read_constructor,
                          size_t total_number_of_reads) :
                config_(config),
                clone_set_(clone_set),
                clone_by_read_constructor_(clone_by_read_constructor),
                total_number_of_reads_(total_number_of_reads) {

            for(int i = 0; i < config.run_params.num_threads; i++)
                thread_tree_storages_.push_back(EvolutionaryTreeStorage());
        }

        EvolutionaryTreeStorage ConstructClonalTrees();
    };
}