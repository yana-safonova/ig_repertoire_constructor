#pragma once

#include "../graph_utils/sparse_graph.hpp"
#include "../graph_utils/decomposition.hpp"

namespace dense_subgraph_finder {

    class DecompositionStatsCalculator {
        DecompositionPtr decomposition_;
        SparseGraphPtr hamming_graph_;

        // auxiliary stats
        std::vector <size_t> num_edges_in_class_;
        std::vector<double> class_edge_fillin_;

        // short stats
        double max_edge_fillin_;
        double fillin_of_max_class_;
        double average_fillin_;
        size_t num_trivial_classes_;

        size_t GetClassOfVertices(size_t v1, size_t v2);

        void Initialize();

        double ComputeClassFillin(size_t class_id);

        void ComputeShortStats();

    public:
        DecompositionStatsCalculator(DecompositionPtr decomposition,
                                     SparseGraphPtr hamming_graph) :
                decomposition_(decomposition),
                hamming_graph_(hamming_graph),
                max_edge_fillin_(0),
                fillin_of_max_class_(0),
                average_fillin_(0),
                num_trivial_classes_(0) {
            Initialize();
        }

        void WriteAllStats(ostream &out);

        void WriteShortStats(ostream &out);

    private:
        DECL_LOGGER("DecompositionStatsCalculator");
    };

}