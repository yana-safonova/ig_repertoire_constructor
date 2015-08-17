#pragma once

#include "../graph_utils/sparse_graph.hpp"
#include "../graph_utils/decomposition.hpp"

namespace dense_subgraph_finder {

    class DecompositionStatsCalculator {
        DecompositionPtr decomposition_;
        SparseGraphPtr hamming_graph_;

        std::vector <size_t> num_edges_in_class;

        size_t GetClassOfVertices(size_t v1, size_t v2);

        void Initialize();

        double ComputeClassFillin(size_t class_id);

        size_t ComputeClassSize(size_t class_id);

    public:
        DecompositionStatsCalculator(DecompositionPtr decomposition,
                                     SparseGraphPtr hamming_graph) :
                decomposition_(decomposition),
                hamming_graph_(hamming_graph) { }

        void WriteStatsInFile(string filename);

    private:
        DECL_LOGGER("DecompositionStatsCalculator");
    };

}