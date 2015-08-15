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

DecompositionPtr MetisDenseSubgraphConstructor::AddIsolatedVertices(SparseGraphPtr hamming_graph_ptr,
                                                                    DecompositionPtr primary_decomposition_ptr) {
    size_t num_isolated_vertices = 0;
    for(size_t i = 0; i < hamming_graph_ptr->N(); i++)
        if(hamming_graph_ptr->VertexIsIsolated(i)) {
            primary_decomposition_ptr->SetClass(i, primary_decomposition_ptr->LastClassId() + 1);
            num_isolated_vertices++;
        }
    TRACE(num_isolated_vertices << " isolated vertices were added")
    return primary_decomposition_ptr;
}

DecompositionPtr MetisDenseSubgraphConstructor::ExpandDecomposition(SparseGraphPtr hamming_graph_ptr,
                                       GraphCollapsedStructurePtr collapsed_struct_ptr,
                                       DecompositionPtr prefinal_decomposition_ptr) {
    DecompositionPtr final_decomposition_ptr1(new Decomposition(hamming_graph_ptr->N()));
    for(size_t i = 0; i < hamming_graph_ptr->N(); i++) {
        size_t main_vertex_id = collapsed_struct_ptr->GetMainVertexIndex(i);
        size_t new_main_index = collapsed_struct_ptr->NewIndexOfOldVertex(main_vertex_id);
        size_t class_id = prefinal_decomposition_ptr->GetVertexClass(new_main_index);
        final_decomposition_ptr1->SetClass(i, class_id);
    }
    DecompositionPtr final_decomposition_ptr2(new Decomposition(hamming_graph_ptr->N()));
    size_t cur_class_id = 0;
    for(size_t i = 0; i < final_decomposition_ptr1->Size(); i++)
        if(final_decomposition_ptr1->ClassSize(i) != 0) {
            auto cur_class = final_decomposition_ptr1->GetClass(i);
            for(auto it = cur_class.begin(); it != cur_class.end(); it++)
                final_decomposition_ptr2->SetClass(*it, cur_class_id);
            cur_class_id++;
        }
    return final_decomposition_ptr2;
}

DecompositionPtr MetisDenseSubgraphConstructor::CreateDecomposition(SparseGraphPtr hamming_graph_ptr,
                                     GraphCollapsedStructurePtr collapsed_struct_ptr) {
    TRACE("== Computation of permutation using METIS");
    TRACE("Input graph contains " << hamming_graph_ptr->N() << " vertices & " << hamming_graph_ptr->NZ() << " edges");
    if(hamming_graph_ptr->NZ() == 0) {
        TRACE("Graph does not contain edges. Trivial decomposition was created");
        DecompositionPtr trivial_decomposition(new Decomposition(hamming_graph_ptr->N()));
        trivial_decomposition = AddIsolatedVertices(hamming_graph_ptr, trivial_decomposition);
        trivial_decomposition->SaveTo(decomposition_filename_);
        return trivial_decomposition;
    }
    PermutationPtr permutation_ptr = CreatePermutation(hamming_graph_ptr, collapsed_struct_ptr);
    TRACE("Computation of primary dense subgraph decomposition starts");
    DecompositionPtr primary_decomposition_ptr = CreatePrimaryDecomposition(hamming_graph_ptr,
                                                                            collapsed_struct_ptr, permutation_ptr);
    TRACE("Primary decomposition contains " << primary_decomposition_ptr->Size() << " subgraphs");
    TRACE("Improvement of the primary decomposition starts");
    DecompositionPtr dense_sgraph_decomposition = ImprovePrimaryDecomposition(hamming_graph_ptr,
                                                                              collapsed_struct_ptr,
                                                                              primary_decomposition_ptr);
    TRACE("Improved decomposition contains " << dense_sgraph_decomposition->Size() << " subgraphs");
    TRACE("Expansion of the constructed decomposition and addition of isolated vertices");
    dense_sgraph_decomposition = ExpandDecomposition(hamming_graph_ptr,
                                                       collapsed_struct_ptr,
                                                       dense_sgraph_decomposition);
    TRACE("Final decomposition contains " << dense_sgraph_decomposition->Size() << " subgraphs");
    dense_sgraph_decomposition->SaveTo(decomposition_filename_);
    return dense_sgraph_decomposition;
}