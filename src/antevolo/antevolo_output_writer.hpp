#pragma once

#include "antevolo_config.hpp"
#include "evolutionary_tree_annotation/annotated_tree_storage.hpp"
#include "evolutionary_graph_utils/evolutionary_tree.hpp"

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

        void WriteTreeInFile(std::string output_dir, const EvolutionaryTree& tree) const;

        void WriteTreeVerticesInFile(std::string output_dir, const EvolutionaryTree& tree) const;


    private:

        void WriteEdge(const EvolutionaryEdgePtr& edge, std::ofstream& out) const;

    };
}