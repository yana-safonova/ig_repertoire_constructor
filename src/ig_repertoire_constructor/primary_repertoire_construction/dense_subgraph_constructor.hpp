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
	double primary_edge_fillin_threshold_;
	double dense_subgraph_joining_threshold_;

	// output struct
	HG_DecompositionPtr dense_subgraph_decomposition_ptr_;

	PermutationPtr CreatePermutation(CRS_HammingGraph_Ptr hamming_graph_ptr,
			HG_CollapsedStructs_Ptr collapsed_struct_ptr, size_t graph_id) {
		return MetisPermutationConstructor(hamming_graph_ptr, collapsed_struct_ptr, graph_id).CreatePermutation();
	}

	HG_DecompositionPtr CreatePrimaryDecomposition(CRS_HammingGraph_Ptr hamming_graph_ptr,
			HG_CollapsedStructs_Ptr collapsed_struct_ptr,
			PermutationPtr permutation_ptr) {
        SimpleDecompositionConstructor simple_constructor(hamming_graph_ptr,
        		permutation_ptr, collapsed_struct_ptr, primary_edge_fillin_threshold_);
        return simple_constructor.CreateDecomposition();

        // todo: introduce config param and compute stats if it is enabled
        //DecompositionStatsCalculator calculator(primary_decomposition_ptr_, hamming_graph_ptr_, collapsed_struct_);
        //calculator.WriteStatsInFile("primary_decomposition_stats.txt");
	}

	HG_DecompositionPtr ImprovePrimaryDecomposition(CRS_HammingGraph_Ptr hamming_graph_ptr,
			HG_CollapsedStructs_Ptr collapsed_struct_ptr,
			HG_DecompositionPtr primary_decomposition_ptr) {
        GreedyJoiningDecomposition decomposition_improver(hamming_graph_ptr, collapsed_struct_ptr,
                primary_decomposition_ptr, dense_subgraph_joining_threshold_);
        return decomposition_improver.ConstructDecomposition();

        // todo: introduce config param and compute stats if it is enabled
        //DecompositionStatsCalculator calculator(dense_subgraph_decomposition_ptr_,
        //		hamming_graph_ptr_, collapsed_struct_ptr_);
        //stringstream ss;
        //ss << "secondary_decomposition_stats_" << dense_subgraph_joining_threshold_ << ".txt";
        //calculator.WriteStatsInFile(ss.str());
	}

public:
	MetisDenseSubgraphConstructor(double primary_edge_fillin_threshold,
			double dense_subgraph_joining_threshold) :
			primary_edge_fillin_threshold_(primary_edge_fillin_threshold),
			dense_subgraph_joining_threshold_(dense_subgraph_joining_threshold) { }

	HG_DecompositionPtr CreateDecomposition(CRS_HammingGraph_Ptr hamming_graph_ptr,
			HG_CollapsedStructs_Ptr collapsed_struct_ptr, size_t graph_id) {
		TRACE("Computation of permutation using METIS");
		PermutationPtr permutation_ptr = CreatePermutation(hamming_graph_ptr, collapsed_struct_ptr, graph_id);
		TRACE("Computation of primary dense subgraph decomposition");
		HG_DecompositionPtr primary_decomposition_ptr = CreatePrimaryDecomposition(hamming_graph_ptr,
				collapsed_struct_ptr, permutation_ptr);
		TRACE("Improvement of the primary decomposition");
		return ImprovePrimaryDecomposition(hamming_graph_ptr, collapsed_struct_ptr,
				primary_decomposition_ptr);
	}

private:
	DECL_LOGGER("MetisDenseSubgraphConstructor");
};

}
