#pragma once

#include "antevolo_config.hpp"
#include "evolutionary_tree_annotation/annotated_tree_storage.hpp"

namespace antevolo {
    class AntEvoloOutputWriter {
        const AntEvoloConfig::OutputParams &output_params_;
        const AnnotatedTreeStorage &annotated_storage_;

    public:
        AntEvoloOutputWriter(const AntEvoloConfig::OutputParams &output_params,
                             const AnnotatedTreeStorage &annotated_storage) :
                output_params_(output_params),
                annotated_storage_(annotated_storage) { }

        void OutputTreeStats() const;

        void OutputSHMForTrees() const;
    };
}