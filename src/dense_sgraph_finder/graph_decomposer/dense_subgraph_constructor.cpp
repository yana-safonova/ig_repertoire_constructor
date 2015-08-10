#include "dense_subgraph_constructor.hpp"

using namespace dense_subgraph_finder;

PermutationPtr MetisDenseSubgraphConstructor::CreatePermutation(SparseGraphPtr hamming_graph_ptr,
                                 GraphCollapsedStructurePtr collapsed_struct_ptr) {
    return MetisPermutationConstructor(hamming_graph_ptr,
                                       collapsed_struct_ptr,
                                       metis_params_,
                                       graph_filename_).CreatePermutation();
}

DecompositionPtr MetisDenseSubgraphConstructor::CreatePrimaryDecomposition(SparseGraphPtr hamming_graph_ptr,
                                            GraphCollapsedStructurePtr collapsed_struct_ptr,
                                            PermutationPtr permutation_ptr) {
    SimpleDecompositionConstructor simple_constructor(hamming_graph_ptr,
                                                      permutation_ptr,
                                                      collapsed_struct_ptr,
                                                      dsf_params_.class_joining_edge_threshold);
    return simple_constructor.CreateDecomposition();
}

DecompositionPtr MetisDenseSubgraphConstructor::ImprovePrimaryDecomposition(SparseGraphPtr hamming_graph_ptr,
                                             GraphCollapsedStructurePtr collapsed_struct_ptr,
                                             DecompositionPtr primary_decomposition_ptr) {
    GreedyJoiningDecomposition decomposition_improver(hamming_graph_ptr, collapsed_struct_ptr,
                                                      primary_decomposition_ptr, dsf_params_.edge_perc_threshold);
    return decomposition_improver.ConstructDecomposition();
}

DecompositionPtr MetisDenseSubgraphConstructor::CreateDecomposition(SparseGraphPtr hamming_graph_ptr,
                                     GraphCollapsedStructurePtr collapsed_struct_ptr) {
    INFO("== Computation of permutation using METIS");
    PermutationPtr permutation_ptr = CreatePermutation(hamming_graph_ptr, collapsed_struct_ptr);
    INFO("== Computation of primary dense subgraph decomposition");
    DecompositionPtr primary_decomposition_ptr = CreatePrimaryDecomposition(hamming_graph_ptr,
                                                                            collapsed_struct_ptr, permutation_ptr);
    INFO("Primary decomposition contains " << primary_decomposition_ptr->Size() << " subgraphs");
    INFO("== Improvement of the primary decomposition");
    DecompositionPtr dense_sgraph_decomposition = ImprovePrimaryDecomposition(hamming_graph_ptr,
                                                                              collapsed_struct_ptr,
                                                                              primary_decomposition_ptr);
    INFO("Final decomposition contains " << dense_sgraph_decomposition->Size() << " subgraphs");
    return dense_sgraph_decomposition;
}