#pragma once

#include <../graph_utils/sparse_graph.hpp>

namespace cdr_labeler {
    class DiversityAnalyser {
        SparseGraphPtr cdr_graph_;

        void InitializeGraph(std::string compressed_cdr3_fasta);

        size_t ComputeD50(const std::vector<SparseGraphPtr> connected_components) const;

    public:
        DiversityAnalyser(std::string compressed_cdr3_fasta) {
            InitializeGraph(compressed_cdr3_fasta);
        }
    };
}