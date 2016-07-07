#pragma once

#include <annotation_utils/annotated_clone_set.hpp>
#include "antevolo_config.hpp"

namespace antevolo {
    class NaiveAntEvoloProcessing {
        const AntEvoloConfig& config_;
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;

        std::string GetTreeOutputFname(std::string output_dir, size_t index1, size_t index2, size_t tree_size);

    public:
        NaiveAntEvoloProcessing(const AntEvoloConfig& config,
                                const annotation_utils::CDRAnnotatedCloneSet &clone_set) :
                config_(config),
                clone_set_(clone_set) { }

        void ConstructClonalTrees();
    };
}