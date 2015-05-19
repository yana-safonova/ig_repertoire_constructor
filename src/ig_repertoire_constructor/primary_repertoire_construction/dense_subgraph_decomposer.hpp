#pragma once

#include "hamming_graph_clusterization/crs_matrix.hpp"
#include "hamming_graph_clusterization/hg_decomposition.hpp"
#include "../spliced_read.hpp"

namespace ig_repertoire_constructor {

class IterativeDenseSubgraphDecomposer {
	size_t min_recessive_abs_size_;
	double min_recessive_rel_size_;

public:
	IterativeDenseSubgraphDecomposer(size_t min_recessive_abs_size, double min_recessive_rel_size) :
		min_recessive_abs_size_(min_recessive_abs_size),
		min_recessive_rel_size_(min_recessive_rel_size) { }

	HG_DecompositionPtr CreateDecomposition(CRS_HammingGraph_Ptr hamming_graph_ptr,
			HG_CollapsedStructs_Ptr collapsed_struct_ptr,
			HG_DecompositionPtr dense_subgraph_decomposition,
			SplicedReadGroup read_group) {
		// todo: implement me!
		return HG_DecompositionPtr(new HG_Decomposition(hamming_graph_ptr->N()));
	}

private:
	DECL_LOGGER("IterativeDenseSubgraphDecomposer");
};

}
