#pragma once

#include "hamming_graph_clusterization/crs_matrix.hpp"
#include "hamming_graph_clusterization/permutation.hpp"
#include "hamming_graph_clusterization/hg_decomposition.hpp"

#include "hamming_graph_clusterization/simple_decomposition_constructor.hpp"
#include "hamming_graph_clusterization/greedy_joining_decomposition_constructor.hpp"
#include "hamming_graph_clusterization/decomposition_stats_calculator.hpp"
#include "metis_permutation_constructor.hpp"

namespace ig_repertoire_constructor {

class MetisDenseSubgraphConstructor {
	// config params
	const ig_config::hg_clusterization_params& params_;

	// output struct
	HG_DecompositionPtr dense_subgraph_decomposition_ptr_;

	PermutationPtr CreatePermutation(CRS_HammingGraph_Ptr hamming_graph_ptr,
			HG_CollapsedStructs_Ptr collapsed_struct_ptr, size_t graph_id) {
		return MetisPermutationConstructor(hamming_graph_ptr, collapsed_struct_ptr, graph_id, params_.hgc_io_params).CreatePermutation();
	}

	HG_DecompositionPtr CreatePrimaryDecomposition(CRS_HammingGraph_Ptr hamming_graph_ptr,
			HG_CollapsedStructs_Ptr collapsed_struct_ptr,
			PermutationPtr permutation_ptr) {
        SimpleDecompositionConstructor simple_constructor(hamming_graph_ptr,
        		permutation_ptr, collapsed_struct_ptr, params_.class_joining_edge_threshold);
        return simple_constructor.CreateDecomposition();
	}

	HG_DecompositionPtr ImprovePrimaryDecomposition(CRS_HammingGraph_Ptr hamming_graph_ptr,
			HG_CollapsedStructs_Ptr collapsed_struct_ptr,
			HG_DecompositionPtr primary_decomposition_ptr) {
        GreedyJoiningDecomposition decomposition_improver(hamming_graph_ptr, collapsed_struct_ptr,
                primary_decomposition_ptr, params_.edge_perc_threshold);
        return decomposition_improver.ConstructDecomposition();
	}

public:
	MetisDenseSubgraphConstructor(const ig_config::hg_clusterization_params& params) :
		params_(params) { }

	HG_DecompositionPtr CreateDecomposition(CRS_HammingGraph_Ptr hamming_graph_ptr,
			HG_CollapsedStructs_Ptr collapsed_struct_ptr, size_t graph_id) {
		TRACE("Computation of permutation using METIS");
		PermutationPtr permutation_ptr = CreatePermutation(hamming_graph_ptr, collapsed_struct_ptr, graph_id);
		TRACE("Computation of primary dense subgraph decomposition");
		HG_DecompositionPtr primary_decomposition_ptr = CreatePrimaryDecomposition(hamming_graph_ptr,
				collapsed_struct_ptr, permutation_ptr);
		TRACE("Improvement of the primary decomposition");
		HG_DecompositionPtr dense_sgraph_decomposition = ImprovePrimaryDecomposition(hamming_graph_ptr, collapsed_struct_ptr,
				primary_decomposition_ptr);
		return dense_sgraph_decomposition;
	}

private:
	DECL_LOGGER("MetisDenseSubgraphConstructor");
};

}
