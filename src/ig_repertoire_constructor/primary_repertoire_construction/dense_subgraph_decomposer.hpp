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

	// aux params;
	size_t best_column_;
	ClassIdRange new_classes_range_;

	typedef map<char, size_t> ColumnMap;

	typedef map<char, vector<size_t> > ColumnIndexMap;

	ColumnIndexMap ComputeColumnIndexMap(size_t pos) {
		ColumnIndexMap cimap;
		for(auto it = decomposition_class_.begin(); it != decomposition_class_.end(); it++) {
			size_t new_index = *it;
			size_t old_index = collapsed_struct_->OldIndexOfNewVertex(new_index);
			char nucleotide = nucl(read_group_[old_index][pos]);
			if(cimap.find(nucleotide) == cimap.end())
				cimap[nucleotide] = vector<size_t>();
			cimap[nucleotide].push_back(new_index);
		}
		return cimap;
	}

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

	size_t MapSize(ColumnMap cmap) {
		size_t map_size = 0;
		for(auto it = cmap.begin(); it != cmap.end(); it++)
			map_size += it->second;
		return map_size;
	}

	bool ColumnIsGood(ColumnMap cmap) {
		size_t recessive_size = GetColumnRecessiveSize(cmap);
		size_t map_size = MapSize(cmap);
		return (recessive_size >= min_recessive_abs_size_) and
				(double(recessive_size) / double(map_size) >= min_recessive_rel_size_);
	}

	size_t UpdateBestColumn(size_t column, ColumnMap cmap, size_t best_recessive_size) {
		size_t recessive_size = GetColumnRecessiveSize(cmap);
		if(recessive_size > best_recessive_size) {
			best_column_ = column;
			return recessive_size;
		}
		return best_recessive_size;
	}

	void FindBestAlignmentColumn() {
		size_t spliced_read_length = read_group_[0].size();
		size_t best_recessive_size = 0;
		size_t num_nt_columns = 0;
		for(size_t i = 0; i < spliced_read_length; i++) {
			auto column_map = ComputeColumnMap(i);
			if(ColumnIsTrivial(column_map))
				continue;
			if(!ColumnIsGood(column_map))
				continue;
			num_nt_columns++;
			TRACE("Column at position " << i << " is good. Map: ");
			for(auto it = column_map.begin(); it != column_map.end(); it++)
				TRACE(it->first << " " << it->second);
			best_recessive_size = UpdateBestColumn(i, column_map, best_recessive_size);
		}
		if(num_nt_columns != 0) {
			TRACE("Best column: " << best_column_ << " with recessive size " <<
				best_recessive_size);
		}
		else {
			TRACE("Not trivial columns were not found");
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
		min_recessive_rel_size_(min_recessive_rel_size),
		best_column_(size_t(-1)) {
		FindBestAlignmentColumn();
	}

	bool ClassCanBeDecomposed() {
		return best_column_ != size_t(-1);
	}

	HG_DecompositionPtr DecomposeClass() {
		assert(best_column_ != size_t(-1));
		ColumnIndexMap cimap = ComputeColumnIndexMap(best_column_);
		size_t new_class_id = output_decomposition_->NextClassId();
		new_classes_range_.first = new_class_id;
		for(auto it = cimap.begin(); it != cimap.end(); it++) {
			auto current_subclass = it->second;
			for(auto elem = current_subclass.begin(); elem != current_subclass.end(); elem++) {
				output_decomposition_->SetClass(*elem, new_class_id);
			}
			TRACE(current_subclass.size() << " elements were added into class #" << new_class_id);
			new_class_id++;
		}
		new_classes_range_.second = new_class_id;
		TRACE("New classes range: " << new_classes_range_.first << " - " << new_classes_range_.second);
		TRACE("# classes in new decomposition: " << output_decomposition_->Size());
		return output_decomposition_;
	}

	ClassIdRange NewClassesIds() {
		return new_classes_range_;
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
	set<size_t> classes_for_analysis_;
	HG_DecompositionPtr current_decomposition_;
	HG_DecompositionPtr new_decomposition_;

	void ClearAuxiliaryStructs() {
		classes_for_analysis_.clear();
	}

//	bool IsDecompositionClassTrivial(const set<size_t> &cur_class) {
//		return cur_class.size() < min_dense_sgraph_size_;
//	}

	bool DecompositionClassIsTrivial(HG_DecompositionPtr decomposition, size_t class_id) {
		return decomposition->RealSizeOfClass(class_id, collapsed_struct_ptr_) < min_dense_sgraph_size_;
	}

	void AddClassIntoNewDecomposition(const DecompositionClass &current_class) {
		size_t new_class_id = new_decomposition_->NextClassId();
		TRACE("Class of size " << current_class.size() << " will be written with id " << new_class_id);
		for(auto it = current_class.cbegin(); it != current_class.cend(); it++)
			new_decomposition_->SetClass(*it, new_class_id);
	}

	void DefineClassesForAnalysis() {
		ClearAuxiliaryStructs();
		for(size_t i = 0; i < current_decomposition_->Size(); i++) {
            if(!DecompositionClassIsTrivial(current_decomposition_, i))
            	classes_for_analysis_.insert(i);
		}
		TRACE(classes_for_analysis_.size() << " classes for analysis were computed");
	}

	void CopyTrivialClasses() {
		for(size_t i = 0; i < current_decomposition_->Size(); i++) {
            auto current_class = current_decomposition_->GetClass(i);
            if(classes_for_analysis_.find(i) == classes_for_analysis_.end())
            	AddClassIntoNewDecomposition(current_class);
		}
	}

	void InitializeNewDecomposition() {
		new_decomposition_ = HG_DecompositionPtr(
				new HG_Decomposition(current_decomposition_->VertexNumber()));
	}

	bool ContinueDecomposition() {
		return classes_for_analysis_.size() > 0;
	}

	ClassIdRange DecomposeClass(const DecompositionClass &current_class) {
		UniformClassDecomposer class_decomposer(current_class,
				new_decomposition_, read_group_, hamming_graph_ptr_, collapsed_struct_ptr_,
				min_recessive_abs_size_, min_recessive_rel_size_);
		if(!class_decomposer.ClassCanBeDecomposed()) {
			TRACE("Class can not be decomposed");
			AddClassIntoNewDecomposition(current_class);
			return ClassIdRange(0, 0);
		}
		new_decomposition_ = class_decomposer.DecomposeClass();
		return class_decomposer.NewClassesIds();
	}

	void TraceNewClasses(const set<size_t> &new_classes_for_analysis) {
		stringstream ss;
		for(auto it = new_classes_for_analysis.cbegin(); it != new_classes_for_analysis.cend(); it++)
			ss << *it << " ";
		TRACE("New classes for analysis: " << ss.str());
	}

	void CheckNewClasses(ClassIdRange new_classes_range, set<size_t> &new_classes_for_analysis) {
		TRACE("Checking new classes");
		for(size_t i = new_classes_range.first; i < new_classes_range.second; i++) {
			auto new_class = new_decomposition_->GetClass(i);
			TRACE("Checking new class #" << i << " of size " << new_class.size());
			if(!DecompositionClassIsTrivial(new_decomposition_, i))
				new_classes_for_analysis.insert(i);
		}

		if(new_classes_for_analysis.size() == 0) {
			TRACE("No new classes for analysis");
			return;
		}
		TraceNewClasses(new_classes_for_analysis);
	}

	void CreateNewDecomposition() {
		TRACE("New iteration of decomposition");
		set<size_t> new_classes_for_analysis;
		TRACE(classes_for_analysis_.size() << " class(es) will be analysed");
		for(auto it = classes_for_analysis_.begin(); it!= classes_for_analysis_.end(); it++) {
			TRACE("== Analysis of old class #" << *it << ", size: " << current_decomposition_->ClassSize(*it));
			auto newly_constructed_classes = DecomposeClass(current_decomposition_->GetClass(*it));
			CheckNewClasses(newly_constructed_classes, new_classes_for_analysis);
		}
		classes_for_analysis_ = new_classes_for_analysis;
	}

	void PrepareNewIteration() {
		current_decomposition_ = new_decomposition_;
		InitializeNewDecomposition();
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
	}

	HG_DecompositionPtr CreateDecomposition() {
		current_decomposition_ = dense_subgraph_decomposition_;
		TRACE("Iterative decomposer of dense subgraph starts");
		size_t num_iteration = 0;
		DefineClassesForAnalysis();
		InitializeNewDecomposition();
		while(ContinueDecomposition()) {
			CopyTrivialClasses();
			CreateNewDecomposition();
			num_iteration++;
			TRACE("----------");
			TRACE("Iteration #" << num_iteration << " finished");
			TRACE("New decomposition contains " << new_decomposition_->Size() << " classes");
			PrepareNewIteration();
		}
		TRACE("Final decomposition consists of " << current_decomposition_->Size() << " classes");
		TRACE("Iterative decomposer of dense subgraph ends");
		return current_decomposition_;
	}

private:
	DECL_LOGGER("IterativeDenseSubgraphDecomposer");
};

}
