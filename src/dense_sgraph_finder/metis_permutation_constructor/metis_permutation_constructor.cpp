#include "metis_permutation_constructor.hpp"

namespace dense_subgraph_finder {

std::string MetisPermutationConstructor::GetMETISGraphFilename() {
    std::stringstream ss;
    ss << path::append_path(params_.hg_output_dir, "hamming_graph");
    ss << "_" << graph_id_ << "_size_" << collapsed_struct_ptr_->NumberNewVertices() << ".graph";
    return ss.str();
}

void MetisPermutationConstructor::WriteHammingGraphInMETISFormat(std::string graph_fname) {
    std::ofstream output_fhandler(graph_fname.c_str());
    output_fhandler << collapsed_struct_ptr_->NumberNewVertices() << "\t" << GetNumEdgesInCollapsedGraph() << endl;
    for (size_t i = 0; i < hamming_graph_ptr_->N(); i++) {
        if (!collapsed_struct_ptr_->VertexIsMain(i))
            continue;
        for (size_t j = hamming_graph_ptr_->RowIndexT()[i]; j < hamming_graph_ptr_->RowIndexT()[i + 1]; j++) {
            size_t v = hamming_graph_ptr_->ColT()[j];
            if (!collapsed_struct_ptr_->VertexIsMain(v))
                continue;
            output_fhandler << collapsed_struct_ptr_->NewIndexOfOldVertex(v) + 1 << "\t";
        }
        for (size_t j = hamming_graph_ptr_->RowIndex()[i]; j < hamming_graph_ptr_->RowIndex()[i + 1]; j++) {
            size_t v = hamming_graph_ptr_->Col()[j];
            if (!collapsed_struct_ptr_->VertexIsMain(v))
                continue;
            output_fhandler << collapsed_struct_ptr_->NewIndexOfOldVertex(v) + 1 << "\t";
        }
        output_fhandler << std::endl;
    }
    output_fhandler.close();
}

std::string MetisPermutationConstructor::RunMETIS(std::string graph_fname) {
    std::string command_line = params_.run_metis + " " + graph_fname + " > " + params_.trash_output;
    int err_code = system(command_line.c_str());
    TRACE("Error code: " << err_code);
    return graph_fname + ".iperm";
}

PermutationPtr MetisPermutationConstructor::ReadPermutation(std::string permutation_fname) {
    PermutationPtr perm = PermutationPtr(new Permutation(collapsed_struct_ptr_->NumCollapsedVertices()));
    perm->ReadFromFile(permutation_fname);
    return perm;
}

PermutationPtr MetisPermutationConstructor::CreatePermutation() {
    std::string metis_graph_fname = GetMETISGraphFilename();
    WriteHammingGraphInMETISFormat(metis_graph_fname);
    std::string permutation_fname = RunMETIS(metis_graph_fname);
    return ReadPermutation(permutation_fname);
}

}