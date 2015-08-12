#pragma once

#include <dsf_config.hpp>
#include "graph_utils/permutation.hpp"
#include "graph_utils/sparse_graph.hpp"
#include "graph_utils/graph_collapsed_structure.hpp"

namespace dense_subgraph_finder {

class MetisPermutationConstructor {
	SparseGraphPtr hamming_graph_ptr_;
	GraphCollapsedStructurePtr collapsed_struct_ptr_;
	const dsf_config::metis_io_params &metis_io_params_;
	const dsf_config::io_params &io_params_;

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
			const dsf_config::metis_io_params &metis_io_params,
            const dsf_config::io_params &io_params) :
		hamming_graph_ptr_(hamming_graph_ptr),
		collapsed_struct_ptr_(collapsed_struct_ptr),
		metis_io_params_(metis_io_params),
		io_params_(io_params) { }

	PermutationPtr CreatePermutation();
};

}
