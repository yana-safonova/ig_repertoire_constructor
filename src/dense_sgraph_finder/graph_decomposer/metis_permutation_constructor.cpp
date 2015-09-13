#include "metis_permutation_constructor.hpp"

using namespace dense_subgraph_finder;

// todo: remove this function
std::string MetisPermutationConstructor::GetMETISGraphFilename() {
    return graph_filename_;
}

void MetisPermutationConstructor::WriteHammingGraphInMETISFormat(std::string graph_fname) {
    std::ofstream output_fhandler(graph_fname.c_str());
    output_fhandler << graph_ptr_->N() << "\t" << graph_ptr_->NZ() << endl;
    for (size_t i = 0; i < graph_ptr_->N(); i++) {
        for (size_t j = graph_ptr_->RowIndexT()[i]; j < graph_ptr_->RowIndexT()[i + 1]; j++) {
            size_t v = graph_ptr_->ColT()[j];
            output_fhandler << v + 1 << "\t";
        }
        for (size_t j = graph_ptr_->RowIndex()[i]; j < graph_ptr_->RowIndex()[i + 1]; j++) {
            size_t v = graph_ptr_->Col()[j];
            output_fhandler << v + 1 << "\t";
        }
        output_fhandler << std::endl;
    }
    output_fhandler.close();
}

std::string MetisPermutationConstructor::RunMETIS(std::string graph_fname) {
    std::string command_line = metis_io_params_.run_metis + " " + graph_fname + " > " + metis_io_params_.trash_output;
    int err_code = system(command_line.c_str());
    string permutation_fname = graph_filename_ + ".iperm";
    TRACE("Error code: " << err_code);
    return permutation_fname;
}

PermutationPtr MetisPermutationConstructor::ReadPermutation(std::string permutation_fname) {
    PermutationPtr perm = PermutationPtr(new Permutation(graph_ptr_->N()));
    perm->ReadFromFile(permutation_fname);
    return perm;
}

PermutationPtr MetisPermutationConstructor::CreatePermutation() {
    string graph_copy_filename = GetMETISGraphFilename();
    WriteHammingGraphInMETISFormat(graph_copy_filename);
    std::string permutation_fname = RunMETIS(graph_copy_filename);
    TRACE("Permutation was written to " << permutation_fname);
    return ReadPermutation(permutation_fname);
}