#pragma once

#include "hamming_graph_clusterization/crs_matrix.hpp"
#include "hamming_graph_clusterization/hg_decomposition.hpp"
#include "../spliced_read.hpp"

namespace ig_repertoire_constructor {

typedef pair<size_t, size_t> ClassIdRange;

class UniformClassDecomposer {
	// input struct
	const DecompositionClass &decomposition_class_;
	HG_DecompositionPtr output_decomposition_;
	const SplicedReadGroup &read_group_;
	CRS_HammingGraph_Ptr hamming_graph_;
	HG_CollapsedStructs_Ptr collapsed_struct_;

	// input params
	size_t min_recessive_abs_size_;
	double min_recessive_rel_size_;

	typedef map<char, size_t> ColumnMap;

	ColumnMap ComputeColumnMap(size_t pos) {
		ColumnMap result_map;
		for(auto it = decomposition_class_.begin(); it != decomposition_class_.end(); it++) {
			size_t new_index = *it;
			size_t old_index = collapsed_struct_->OldIndexOfNewVertex(new_index);
			size_t multiplicity = collapsed_struct_->MultiplicityOfOldVertex(old_index);
			char nucleotide = nucl(read_group_[old_index][pos]);
			if(result_map.find(nucleotide) == result_map.end())
				result_map[nucleotide] = 0;
			result_map[nucleotide] += multiplicity;
		}
		return result_map;
	}

	size_t GetColumnRecessiveSize(ColumnMap cmap) {
		if(cmap.size() <= 1)
			return 0;
		set<size_t> sizes;
		for(auto it = cmap.begin(); it != cmap.end(); it++)
			sizes.insert(it->second);
		return *(prev(prev(sizes.end())));
	}

	bool ColumnIsTrivial(ColumnMap cmap) {
		if(cmap.size() == 1)
			return true;
		size_t rec_size = GetColumnRecessiveSize(cmap);
		return rec_size <= 1;
	}

	void FindBestAlignmentColumn() {
		size_t spliced_read_length = read_group_[0].size();
		for(size_t i = 0; i < spliced_read_length; i++) {
			auto column_map = ComputeColumnMap(i);
			if(ColumnIsTrivial(column_map))
				continue;
			TRACE("Position: " << i << ", map: ");
			for(auto it = column_map.begin(); it != column_map.end(); it++)
				TRACE(it->first << " " << it->second);
		}
	}

public:
	UniformClassDecomposer(const DecompositionClass &decomposition_class,
			HG_DecompositionPtr output_decomposition,
			const SplicedReadGroup &read_group,
			CRS_HammingGraph_Ptr hamming_graph,
			HG_CollapsedStructs_Ptr collapsed_struct,
			size_t min_recessive_abs_size,
			double min_recessive_rel_size) :
		decomposition_class_(decomposition_class),
		output_decomposition_(output_decomposition),
		read_group_(read_group),
		hamming_graph_(hamming_graph),
		collapsed_struct_(collapsed_struct),
		min_recessive_abs_size_(min_recessive_abs_size),
		min_recessive_rel_size_(min_recessive_rel_size) {
		FindBestAlignmentColumn();
	}

	bool ClassCanBeDecomposed() {
		// todo: imlement me!
		return true;
	}

	HG_DecompositionPtr DecomposeClass() {
		// todo
		return output_decomposition_;
	}

	ClassIdRange NewClassesIds() {
		return ClassIdRange(0, 0);
	}

private:
	DECL_LOGGER("UniformClassDecomposer");
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
		return cur_class.size() < min_dense_sgraph_size_;
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

	ClassIdRange DecomposeClass(const DecompositionClass &current_class, HG_DecompositionPtr new_decomposition) {
		UniformClassDecomposer class_decomposer(current_class,
				new_decomposition, read_group_, hamming_graph_ptr_, collapsed_struct_ptr_,
				min_recessive_abs_size_, min_recessive_rel_size_);
		if(!class_decomposer.ClassCanBeDecomposed()) {
			TRACE("Class can not be decomposed");
			return ClassIdRange(0, 0);
		}
		class_decomposer.DecomposeClass();
		return class_decomposer.NewClassesIds();
	}

	void CheckNewClasses(HG_DecompositionPtr new_decomposition, ClassIdRange new_classes_range,
			vector<size_t> &new_classes_for_analysis) {
		for(size_t i = new_classes_range.first; i < new_classes_range.second; i++) {
			auto new_class = new_decomposition->GetClass(i);
			if(!IsDecompositionClassTrivial(new_class))
				new_classes_for_analysis.push_back(i);
		}
	}

	HG_DecompositionPtr CreateNewDecomposition(HG_DecompositionPtr current_decomposition) {
		TRACE("New iteration of decomposition");
		HG_DecompositionPtr new_decomposition = HG_DecompositionPtr(
				new HG_Decomposition(current_decomposition->VertexNumber()));
		vector<size_t> new_classes_for_analysis;
		TRACE(classes_for_analysis_.size() << " classes will be analysed");
		for(auto it = classes_for_analysis_.begin(); it!= classes_for_analysis_.end(); it++) {
			TRACE("Analysis of class #" << *it << ", size: " << current_decomposition->ClassSize(*it));
			auto newly_constructed_classes = DecomposeClass(
					current_decomposition->GetClass(*it), new_decomposition);
			CheckNewClasses(new_decomposition, newly_constructed_classes, new_classes_for_analysis);
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
		TRACE("Iterative decomposer of dense subgraph starts");
		while(ContinueDecomposition()) {
			new_decomposition = CreateNewDecomposition(current_decomposition);
		}
		TRACE("Final decomposition consists of " << new_decomposition->Size() << " classes");
		TRACE("Iterative decomposer of dense subgraph ends");
		assert(false);
		return new_decomposition;
	}

private:
	DECL_LOGGER("IterativeDenseSubgraphDecomposer");
};

}
