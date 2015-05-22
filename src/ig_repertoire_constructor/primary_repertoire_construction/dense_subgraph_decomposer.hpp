#pragma once

#include "hamming_graph_clusterization/crs_matrix.hpp"
#include "hamming_graph_clusterization/hg_decomposition.hpp"
#include "../spliced_read.hpp"

namespace ig_repertoire_constructor {

class UniformClassDecomposer {
	// input struct
	const DecompositionClass &decomposition_class_;
	HG_DecompositionPtr output_decomposition_;

	// input params
	size_t min_recessive_abs_size_;
	double min_recessive_rel_size_;
	size_t min_dense_sgraph_size_;

public:
	UniformClassDecomposer(const DecompositionClass &decomposition_class,
			HG_DecompositionPtr output_decomposition,
			size_t min_recessive_abs_size,
			double min_recessive_rel_size,
			size_t min_dense_sgraph_size) :
		decomposition_class_(decomposition_class),
		output_decomposition_(output_decomposition),
		min_recessive_abs_size_(min_recessive_abs_size),
		min_recessive_rel_size_(min_recessive_rel_size),
		min_dense_sgraph_size_(min_dense_sgraph_size) { }

	bool ClassCanBeDecomposed() {
		// todo: imlement me!
		return true;
	}

	HG_DecompositionPtr DecomposeClass() {
		// todo
		return output_decomposition_;
	}

	vector<size_t> NewNotTrivialClassesIds() {
		return vector<size_t>();
	}
};

class IterativeDenseSubgraphDecomposer {
	// input params
	size_t min_recessive_abs_size_;
	double min_recessive_rel_size_;
	size_t min_dense_sgraph_size_;

	// input structs
	CRS_HammingGraph_Ptr hamming_graph_ptr_;
	HG_CollapsedStructs_Ptr collapsed_struct_ptr_;
	HG_DecompositionPtr dense_subgraph_decomposition_;
	const SplicedReadGroup &read_group_;

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

	vector<size_t> DecomposeClass(const DecompositionClass &current_class, HG_DecompositionPtr new_decomposition) {
		UniformClassDecomposer class_decomposer(current_class, new_decomposition, min_recessive_abs_size_, min_recessive_rel_size_, min_dense_sgraph_size_);
		if(!class_decomposer.ClassCanBeDecomposed())
			return vector<size_t>();
		class_decomposer.DecomposeClass();
		return class_decomposer.NewNotTrivialClassesIds();
	}

	HG_DecompositionPtr CreateNewDecomposition(HG_DecompositionPtr current_decomposition) {
		HG_DecompositionPtr new_decomposition = HG_DecompositionPtr(new HG_Decomposition(current_decomposition->VertexNumber()));
		vector<size_t> new_classes_for_analysis;
		for(auto it = classes_for_analysis_.begin(); it!= classes_for_analysis_.end(); it++) {
			auto newly_constructed_classes = DecomposeClass(current_decomposition->GetClass(*it), new_decomposition);
			new_classes_for_analysis.insert(new_classes_for_analysis.begin(), newly_constructed_classes.begin(), newly_constructed_classes.end());
		}
		classes_for_analysis_ = new_classes_for_analysis;
		return new_decomposition;
	}

public:
	IterativeDenseSubgraphDecomposer(size_t min_recessive_abs_size,
			double min_recessive_rel_size,
			CRS_HammingGraph_Ptr hamming_graph_ptr,
			HG_CollapsedStructs_Ptr collapsed_struct_ptr,
			HG_DecompositionPtr dense_subgraph_decomposition,
			const SplicedReadGroup &read_group) :
		min_recessive_abs_size_(min_recessive_abs_size),
		min_recessive_rel_size_(min_recessive_rel_size),
		min_dense_sgraph_size_(min_recessive_abs_size * 2 + 1),
		hamming_graph_ptr_(hamming_graph_ptr),
		collapsed_struct_ptr_(collapsed_struct_ptr),
		dense_subgraph_decomposition_(dense_subgraph_decomposition),
		read_group_(read_group) {
		InitializeAuxiliaryStructs(dense_subgraph_decomposition_);
	}

	HG_DecompositionPtr CreateDecomposition() {
		HG_DecompositionPtr current_decomposition = dense_subgraph_decomposition_;
		HG_DecompositionPtr new_decomposition;
		while(ContinueDecomposition()) {
			new_decomposition = CreateNewDecomposition(current_decomposition);
		}
		return new_decomposition;
	}

private:
	DECL_LOGGER("IterativeDenseSubgraphDecomposer");
};

}
