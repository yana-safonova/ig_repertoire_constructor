#pragma once

#include "../graph_utils/sparse_graph.hpp"
#include "../graph_utils/permutation.hpp"
#include "../graph_utils/decomposition.hpp"

#include "simple_decomposition_constructor.hpp"
#include "greedy_joining_decomposition_constructor.hpp"
#include "decomposition_stats_calculator.hpp"
#include "metis_permutation_constructor.hpp"

namespace dense_subgraph_finder {

    class MetisDenseSubgraphConstructor {
        // config params
        const dsf_config::dense_sgraph_finder_params &dsf_params_;
        const dsf_config::metis_io_params &metis_params_;
        string graph_filename_;

        PermutationPtr CreatePermutation(SparseGraphPtr graph_ptr);

        DecompositionPtr CreateDecompositionForSmallGraph(SparseGraphPtr hamming_graph_ptr);

        DecompositionPtr CreatePrimaryDecomposition(SparseGraphPtr hamming_graph_ptr,
                                                    PermutationPtr permutation_ptr);

        DecompositionPtr ImprovePrimaryDecomposition(SparseGraphPtr hamming_graph_ptr,
                                                     DecompositionPtr primary_decomposition_ptr);

        DecompositionPtr AddIsolatedVertices(SparseGraphPtr hamming_graph_ptr,
                                             DecompositionPtr primary_decomposition_ptr);

        DecompositionPtr ExpandDecomposition(SparseGraphPtr hamming_graph_ptr,
                                             DecompositionPtr primary_decomposition_ptr);

    public:
        MetisDenseSubgraphConstructor(const dsf_config::dense_sgraph_finder_params &dsf_params,
                                      const dsf_config::metis_io_params &metis_params,
                                      string graph_filename) :
                dsf_params_(dsf_params),
                metis_params_(metis_params),
                graph_filename_(graph_filename) { }

        DecompositionPtr CreateDecomposition(SparseGraphPtr hamming_graph_ptr);

    private:
        DECL_LOGGER("DenseSubgraphConstructor");
    };

}
