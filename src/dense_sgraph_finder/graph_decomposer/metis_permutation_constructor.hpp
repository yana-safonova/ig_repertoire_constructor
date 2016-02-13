#pragma once

#include "../dsf_config.hpp"
#include "../graph_utils/permutation.hpp"
#include "../graph_utils/sparse_graph.hpp"

namespace dense_subgraph_finder {

class MetisPermutationConstructor {
	SparseGraphPtr graph_ptr_;
	const dsf_config::metis_io_params &metis_io_params_;
	string graph_filename_;

	std::string GetMETISGraphFilename();

	void WriteHammingGraphInMETISFormat(std::string graph_fname);

	std::string RunMETIS(std::string graph_fname);

	PermutationPtr ReadPermutation(std::string permutation_fname);

public:
	PermutationPtr CreatePermutation();

	MetisPermutationConstructor(SparseGraphPtr graph_ptr,
								const dsf_config::metis_io_params &metis_io_params,
								string graph_filename) :
			graph_ptr_(graph_ptr),
			metis_io_params_(metis_io_params),
			graph_filename_(graph_filename) { }
};

}
