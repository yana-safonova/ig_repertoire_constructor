#pragma once

#include <dsf_config.hpp>
#include "graph_utils/permutation.hpp"
#include "graph_utils/sparse_graph.hpp"
#include "graph_utils/graph_collapsed_structure.hpp"

namespace dense_subgraph_finder {

class MetisPermutationConstructor {
	SparseGraphPtr hamming_graph_ptr_;
	GraphCollapsedStructurePtr collapsed_struct_ptr_;
	const dsf_config::metis_io_params &params_;
    string graph_filename;

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
			const dsf_config::metis_io_params &params,
            string graph_filename) :
		hamming_graph_ptr_(hamming_graph_ptr),
		collapsed_struct_ptr_(collapsed_struct_ptr),
		params_(params),
        graph_filename(graph_filename) { }

	PermutationPtr CreatePermutation();
};

}
