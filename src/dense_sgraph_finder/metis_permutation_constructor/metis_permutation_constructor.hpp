#pragma once

#include "hamming_graph_clusterization/permutation.hpp"
#include "graph_utils/crs_matrix.hpp"

namespace ig_repertoire_constructor {

class MetisPermutationConstructor {
	CRS_HammingGraph_Ptr hamming_graph_ptr_;
	HG_CollapsedStructs_Ptr collapsed_struct_ptr_;
	size_t graph_id_;
	const ig_config::hg_clusterization_params::hg_clusterization_io_params &params_;

	std::string GetMETISGraphFilename();

	size_t GetNumEdgesInCollapsedGraph() {
		return collapsed_struct_ptr_->NumberCollapsedEdges(hamming_graph_ptr_);
	}

	void WriteHammingGraphInMETISFormat(std::string graph_fname);

	std::string RunMETIS(std::string graph_fname);

	PermutationPtr ReadPermutation(std::string permutation_fname);

public:
	MetisPermutationConstructor(CRS_HammingGraph_Ptr hamming_graph_ptr,
			HG_CollapsedStructs_Ptr collapsed_struct_ptr,
			size_t graph_id,
			const ig_config::hg_clusterization_params::hg_clusterization_io_params &params) :
		hamming_graph_ptr_(hamming_graph_ptr),
		collapsed_struct_ptr_(collapsed_struct_ptr),
		graph_id_(graph_id),
		params_(params) { }

	PermutationPtr CreatePermutation();
};

}
