#pragma once

#include "hamming_graph_clusterization/permutation.hpp"
#include "hamming_graph_clusterization/crs_matrix.hpp"

namespace ig_repertoire_constructor {

class MetisPermutationConstructor {
	CRS_HammingGraph_Ptr hamming_graph_ptr_;
	HG_CollapsedStructs_Ptr collapsed_struct_ptr_;
	size_t graph_id_;

	// get filename of METIS graph
	string GetMETISGraphFilename() {
		stringstream ss;
		ss << "hamming_graph_" << graph_id_ << ".graph";
		return ss.str();
	}

	size_t GetNumEdgesInCollapsedGraph() {
		// todo: compute me
		return 0;
	}

	// write graph in METIS format
	void WriteHammingGraphInMETISFormat(string graph_fname) {
		ofstream output_fhandler(graph_fname.c_str());
		output_fhandler << hamming_graph_ptr_->N() << "\t" << GetNumEdgesInCollapsedGraph() << endl;
		for(size_t i = 0; i < hamming_graph_ptr_->N(); i++) {
			for(size_t j = hamming_graph_ptr_->RowIndexT()[i]; j < hamming_graph_ptr_->RowIndexT()[i + 1]; j++)
				output_fhandler << hamming_graph_ptr_->ColT()[j] << "\t";
			for(size_t j = hamming_graph_ptr_->RowIndex()[i]; j < hamming_graph_ptr_->RowIndex()[i + 1]; j++)
				output_fhandler << hamming_graph_ptr_->Col()[j] << "\t";
			output_fhandler << endl;
		}
		output_fhandler.close();
	}

	void RunMETIS(string graph_fname) {

	}


public:
	MetisPermutationConstructor(CRS_HammingGraph_Ptr hamming_graph_ptr, HG_CollapsedStructs_Ptr collapsed_struct_ptr, size_t graph_id) :
		hamming_graph_ptr_(hamming_graph_ptr),
		collapsed_struct_ptr_(collapsed_struct_ptr),
		graph_id_(graph_id) { }

	PermutationPtr CreatePermutation() {
		string metis_graph_fname = GetMETISGraphFilename();
		WriteHammingGraphInMETISFormat(metis_graph_fname);
		assert(false);
		return PermutationPtr(new Permutation(collapsed_struct_ptr_->NumCollapsedVertices()));
	}
};

}
