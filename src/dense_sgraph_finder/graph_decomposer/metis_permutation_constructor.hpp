#pragma once

#include <dsf_config.hpp>
#include "graph_utils/permutation.hpp"
#include "graph_utils/sparse_graph.hpp"
#include "graph_utils/graph_collapsed_structure.hpp"

namespace dense_subgraph_finder {

class MetisPermutationConstructor {
	SparseGraphPtr hamming_graph_ptr_;
	GraphCollapsedStructurePtr collapsed_struct_ptr_;
	size_t graph_id_;
	const dsf_config::dense_sgraph_finder_params::hg_clusterization_io_params &params_;

	std::string GetMETISGraphFilename();

	size_t GetNumEdgesInCollapsedGraph() {
		return collapsed_struct_ptr_->NumberCollapsedEdges(hamming_graph_ptr_);
	}

	void WriteHammingGraphInMETISFormat(std::string graph_fname);

	std::string RunMETIS(std::string graph_fname);

	PermutationPtr ReadPermutation(std::string permutation_fname);

public:
	MetisPermutationConstructor(SparseGraphPtr hamming_graph_ptr,
			GraphCollapsedStructurePtr collapsed_struct_ptr,
			size_t graph_id,
			const dsf_config::dense_sgraph_finder_params::hg_clusterization_io_params &params) :
		hamming_graph_ptr_(hamming_graph_ptr),
		collapsed_struct_ptr_(collapsed_struct_ptr),
		graph_id_(graph_id),
		params_(params) { }

	PermutationPtr CreatePermutation();
};

}
