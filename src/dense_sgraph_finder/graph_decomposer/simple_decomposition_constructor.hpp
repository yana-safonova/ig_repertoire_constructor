#pragma once

#include "../graph_utils/sparse_graph.hpp"
#include "../graph_utils/decomposition.hpp"
#include "../graph_utils/permutation.hpp"

namespace dense_subgraph_finder {

    class SimpleDecompositionConstructor {
        // input parameters
        SparseGraphPtr graph_ptr_;
        PermutationPtr permutation_prt_;
        double edge_perc_threshold_;

        //output parameters
        DecompositionPtr decomposition_ptr_;

        void CreateFirstSet();

        double ComputeEdgePercToPreviousSet(size_t vertex);

    public:
        SimpleDecompositionConstructor(SparseGraphPtr graph_ptr,
                                       PermutationPtr permutation_ptr,
                                       double edge_perc_threshold) :
                graph_ptr_(graph_ptr),
                permutation_prt_(permutation_ptr),
                edge_perc_threshold_(edge_perc_threshold),
                decomposition_ptr_(new Decomposition(permutation_ptr->Size())) { }

        DecompositionPtr CreateDecomposition();

    private:
        DECL_LOGGER("SimpleDecompositionConstructor");
    };

}