#include "dense_subgraph_constructor.hpp"

using namespace dense_subgraph_finder;

// if size of input graph is small, we create trivial decomposition:
// all heavy vertices will be located to separate decomposition classes
// light vertiсes will be glued to the first heavy vertex
DecompositionPtr MetisDenseSubgraphConstructor::CreateDecompositionForSmallGraph(SparseGraphPtr hamming_graph_ptr) {
    DecompositionPtr small_graph_decomposition(new Decomposition(hamming_graph_ptr->N()));
    size_t current_class_id = 0;
    for(size_t i = 0; i < hamming_graph_ptr->N(); i++)
        if(hamming_graph_ptr->WeightOfVertex(i) >= dsf_params_.min_supernode_size) {
            small_graph_decomposition->SetClass(i, current_class_id);
            current_class_id += 1;
        }
        else
            small_graph_decomposition->SetClass(i, 0);
    return small_graph_decomposition;
}

PermutationPtr MetisDenseSubgraphConstructor::CreatePermutation(SparseGraphPtr graph_ptr) {
    return MetisPermutationConstructor(graph_ptr,
                                       metis_params_,
                                       graph_filename_).CreatePermutation();
}

DecompositionPtr MetisDenseSubgraphConstructor::CreatePrimaryDecomposition(SparseGraphPtr hamming_graph_ptr,
                                                                           PermutationPtr permutation_ptr) {
    SimpleDecompositionConstructor simple_constructor(hamming_graph_ptr,
                                                      permutation_ptr,
                                                      dsf_params_.primary_edge_fillin,
                                                      dsf_params_.min_supernode_size);
    return simple_constructor.CreateDecomposition();
}

DecompositionPtr MetisDenseSubgraphConstructor::ImprovePrimaryDecomposition(SparseGraphPtr hamming_graph_ptr,
                                                                            DecompositionPtr primary_decomposition_ptr) {
    GreedyJoiningDecomposition decomposition_improver(hamming_graph_ptr,
                                                      primary_decomposition_ptr,
                                                      dsf_params_.min_fillin_threshold,
                                                      dsf_params_.min_supernode_size);
    return decomposition_improver.ConstructDecomposition();
}

DecompositionPtr MetisDenseSubgraphConstructor::AddIsolatedVertices(SparseGraphPtr hamming_graph_ptr,
                                                                    DecompositionPtr primary_decomposition_ptr) {
    size_t num_isolated_vertices = 0;
    for(size_t i = 0; i < hamming_graph_ptr->N(); i++)
        if(hamming_graph_ptr->VertexIsIsolated(i)) {
            primary_decomposition_ptr->SetClass(i, primary_decomposition_ptr->LastClassId() + 1);
            num_isolated_vertices++;
        }
    TRACE(num_isolated_vertices << " isolated vertices were added");
    return primary_decomposition_ptr;
}

DecompositionPtr MetisDenseSubgraphConstructor::ExpandDecomposition(SparseGraphPtr hamming_graph_ptr,
                                                                    DecompositionPtr prefinal_decomposition_ptr) {
    // addition of isolated vertices if they were not processed by METIS
    for(size_t i = 0; i < hamming_graph_ptr->N(); i++) {
        if(!prefinal_decomposition_ptr->VertexClassIsInitialized(i))
            prefinal_decomposition_ptr->SetClass(i, prefinal_decomposition_ptr->NextClassId());
    }
    // assignment consequetive ids for vertices
    DecompositionPtr final_decomposition_ptr(new Decomposition(hamming_graph_ptr->N()));
    size_t cur_class_id = 0;
    for(size_t i = 0; i < prefinal_decomposition_ptr->Size(); i++)
        if(prefinal_decomposition_ptr->ClassSize(i) != 0) {
            auto cur_class = prefinal_decomposition_ptr->GetClass(i);
            for(auto it = cur_class.begin(); it != cur_class.end(); it++)
                final_decomposition_ptr->SetClass(*it, cur_class_id);
            cur_class_id++;
        }
    return final_decomposition_ptr;
}

DecompositionPtr MetisDenseSubgraphConstructor::CreateDecomposition(SparseGraphPtr hamming_graph_ptr) {
    TRACE("== Computation of permutation using METIS");
    TRACE("Input graph contains " << hamming_graph_ptr->N() << " vertices & " << hamming_graph_ptr->NZ() << " edges");
    if(dsf_params_.create_trivial_decomposition or hamming_graph_ptr->N() < dsf_params_.min_graph_size) {
        TRACE("Graph is trivial. Trivial decomposition was created");
        return CreateDecompositionForSmallGraph(hamming_graph_ptr);
    }
    PermutationPtr permutation_ptr = CreatePermutation(hamming_graph_ptr);
    TRACE("Computation of primary dense subgraph decomposition starts");
    DecompositionPtr primary_decomposition_ptr = CreatePrimaryDecomposition(hamming_graph_ptr,
                                                                            permutation_ptr);
    TRACE("Primary decomposition contains " << primary_decomposition_ptr->Size() << " subgraphs");
    TRACE("Improvement of the primary decomposition starts");
    DecompositionPtr dense_sgraph_decomposition = ImprovePrimaryDecomposition(hamming_graph_ptr,
                                                                              primary_decomposition_ptr);
    TRACE("Improved decomposition contains " << dense_sgraph_decomposition->Size() << " subgraphs");
    TRACE("Expansion of the constructed decomposition and addition of isolated vertices");
    dense_sgraph_decomposition = ExpandDecomposition(hamming_graph_ptr,
                                                       dense_sgraph_decomposition);
    TRACE("Final decomposition contains " << dense_sgraph_decomposition->Size() << " subgraphs");
    return dense_sgraph_decomposition;
}