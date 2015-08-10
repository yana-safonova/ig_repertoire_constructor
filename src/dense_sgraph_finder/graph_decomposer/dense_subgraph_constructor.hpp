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
        const dsf_config::dense_sgraph_finder_params &params_;

        // output struct
        DecompositionPtr dense_subgraph_decomposition_ptr_;

        PermutationPtr CreatePermutation(SparseGraphPtr hamming_graph_ptr,
                                         GraphCollapsedStructurePtr collapsed_struct_ptr, size_t graph_id);

        DecompositionPtr CreatePrimaryDecomposition(SparseGraphPtr hamming_graph_ptr,
                                                    GraphCollapsedStructurePtr collapsed_struct_ptr,
                                                    PermutationPtr permutation_ptr);

        DecompositionPtr ImprovePrimaryDecomposition(SparseGraphPtr hamming_graph_ptr,
                                                     GraphCollapsedStructurePtr collapsed_struct_ptr,
                                                     DecompositionPtr primary_decomposition_ptr);

    public:
        MetisDenseSubgraphConstructor(const dsf_config::dense_sgraph_finder_params &params) :
                params_(params) { }

        DecompositionPtr CreateDecomposition(SparseGraphPtr hamming_graph_ptr,
                                             GraphCollapsedStructurePtr collapsed_struct_ptr,
                                             size_t graph_id);

    private:
        DECL_LOGGER("MetisDenseSubgraphConstructor");
    };

}
