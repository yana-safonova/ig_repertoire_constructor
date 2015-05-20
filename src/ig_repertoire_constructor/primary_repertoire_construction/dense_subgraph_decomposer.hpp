#pragma once

#include "hamming_graph_clusterization/crs_matrix.hpp"
#include "hamming_graph_clusterization/hg_decomposition.hpp"
#include "../spliced_read.hpp"

namespace ig_repertoire_constructor {

class IterativeDenseSubgraphDecomposer {
	size_t min_recessive_abs_size_;
	double min_recessive_rel_size_;
	size_t min_dense_sgraph_size_;

	// auxiliary structures
	vector<size_t> classes_for_analysis_;

	void ClearAuxiliaryStructs() {
		classes_for_analysis_.clear();
	}

	bool IsDecompositionClassTrivial(const set<size_t> &cur_class) {
		return cur_class.size() >= min_dense_sgraph_size_;
	}

	void InitializeAuxiliaryStructs(HG_DecompositionPtr primary_decomposition) {
		ClearAuxiliaryStructs();
		for(size_t i = 0; i < primary_decomposition->Size(); i++) {
            auto current_class = primary_decomposition->GetClass(i);
            if(!IsDecompositionClassTrivial(current_class))
            	classes_for_analysis_.push_back(i);
		}
	}

	bool ContinueDecomposition() {
		return classes_for_analysis_.size() > 0;
	}

public:
	IterativeDenseSubgraphDecomposer(size_t min_recessive_abs_size,
			double min_recessive_rel_size,
			size_t min_dense_sgraph_size) :
		min_recessive_abs_size_(min_recessive_abs_size),
		min_recessive_rel_size_(min_recessive_rel_size),
		min_dense_sgraph_size_(min_dense_sgraph_size) { }

	HG_DecompositionPtr CreateDecomposition(CRS_HammingGraph_Ptr hamming_graph_ptr,
			HG_CollapsedStructs_Ptr collapsed_struct_ptr,
			HG_DecompositionPtr dense_subgraph_decomposition,
			SplicedReadGroup read_group) {
		InitializeAuxiliaryStructs(dense_subgraph_decomposition);
		HG_DecompositionPtr current_decomposition_ptr = dense_subgraph_decomposition;
		while(ContinueDecomposition()) {
			for(auto it = classes_for_analysis_.begin(); it != classes_for_analysis_.end(); it++) {
				size_t current_class = *it;
			}
		}
		return current_decomposition_ptr;
	}

private:
	DECL_LOGGER("IterativeDenseSubgraphDecomposer");
};

}
