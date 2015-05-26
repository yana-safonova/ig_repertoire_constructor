#pragma once

#include "hamming_graph_clusterization/permutation.hpp"
#include "hamming_graph_clusterization/crs_matrix.hpp"

namespace ig_repertoire_constructor {

class MetisPermutationConstructor {
	CRS_HammingGraph_Ptr hamming_graph_ptr_;
	HG_CollapsedStructs_Ptr collapsed_struct_ptr_;
	size_t graph_id_;
	const ig_config::hg_clusterization_params::hg_clusterization_io_params &params_;

	// get filename of METIS graph
	string GetMETISGraphFilename() {
		stringstream ss;
		ss << path::append_path(params_.hg_output_dir, "hamming_graph");
		ss << "_" << graph_id_ << "_size_" << collapsed_struct_ptr_->NumberNewVertices() << ".graph";
		return ss.str();
	}

	size_t GetNumEdgesInCollapsedGraph() {
		return collapsed_struct_ptr_->NumberCollapsedEdges(hamming_graph_ptr_);
	}

	// write graph in METIS format
	void WriteHammingGraphInMETISFormat(string graph_fname) {
		ofstream output_fhandler(graph_fname.c_str());
		output_fhandler << collapsed_struct_ptr_->NumberNewVertices() << "\t" << GetNumEdgesInCollapsedGraph() << endl;
		for(size_t i = 0; i < hamming_graph_ptr_->N(); i++) {
			if(!collapsed_struct_ptr_->VertexIsMain(i))
				continue;
			for(size_t j = hamming_graph_ptr_->RowIndexT()[i]; j < hamming_graph_ptr_->RowIndexT()[i + 1]; j++) {
				size_t v = hamming_graph_ptr_->ColT()[j];
				if(!collapsed_struct_ptr_->VertexIsMain(v))
					continue;
				output_fhandler << collapsed_struct_ptr_->NewIndexOfOldVertex(v) + 1 << "\t";
			}
			for(size_t j = hamming_graph_ptr_->RowIndex()[i]; j < hamming_graph_ptr_->RowIndex()[i + 1]; j++) {
				size_t v = hamming_graph_ptr_->Col()[j];
				if(!collapsed_struct_ptr_->VertexIsMain(v))
					continue;
				output_fhandler << collapsed_struct_ptr_->NewIndexOfOldVertex(v) + 1 << "\t";
			}
			output_fhandler << endl;
		}
		output_fhandler.close();
	}

	string RunMETIS(string graph_fname) {
		string command_line = params_.run_metis + " " + graph_fname + " > " + params_.trash_output;
		int err_code = system(command_line.c_str());
		TRACE("Error code: " << err_code);
		return graph_fname + ".iperm";
	}

	PermutationPtr ReadPermutation(string permutation_fname) {
		PermutationPtr perm = PermutationPtr(new Permutation(collapsed_struct_ptr_->NumCollapsedVertices()));
		perm->ReadFromFile(permutation_fname);
		return perm;
	}

public:
	MetisPermutationConstructor(CRS_HammingGraph_Ptr hamming_graph_ptr,
			HG_CollapsedStructs_Ptr collapsed_struct_ptr,
			size_t graph_id,
			const ig_config::hg_clusterization_params::hg_clusterization_io_params &params) :
		hamming_graph_ptr_(hamming_graph_ptr),
		collapsed_struct_ptr_(collapsed_struct_ptr),
		graph_id_(graph_id),
		params_(params) { }

	PermutationPtr CreatePermutation() {
		string metis_graph_fname = GetMETISGraphFilename();
		WriteHammingGraphInMETISFormat(metis_graph_fname);
		string permutation_fname = RunMETIS(metis_graph_fname);
		return ReadPermutation(permutation_fname);
	}
};

}
