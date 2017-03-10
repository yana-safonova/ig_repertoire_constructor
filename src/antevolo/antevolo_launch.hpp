#pragma once

#include "antevolo_config.hpp"
#include "annotation_utils/annotated_clone_set.hpp"

namespace antevolo {
    class AntEvoloLaunch {
        const AntEvoloConfig& config_;

    private:
        void ShmModelPosteriorCalculation(const annotation_utils::AnnotatedCloneSet<annotation_utils::AnnotatedClone>&);

    public:
        AntEvoloLaunch(const AntEvoloConfig& config) : config_(config) { }

        void Launch();
    };
}