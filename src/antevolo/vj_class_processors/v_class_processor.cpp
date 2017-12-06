#include "v_class_processor.hpp"
#include <cdr3_hamming_graph_connected_components_processors/edmonds_cdr3_hg_cc_processor.hpp>


namespace antevolo {

    std::string VClassProcessor::WriteUniqueCDR3InFasta(core::DecompositionClass decomposition_class) {
        std::string output_fname = GetFastaFname(decomposition_class);
        std::ofstream out(output_fname);
        for(auto it = unique_cdr3s_.begin(); it != unique_cdr3s_.end(); it++)
            out << ">" << *it << std::endl << *it << std::endl;
        return output_fname;
    }

    vector<SparseGraphPtr> VClassProcessor::ComputeConnectedComponents(
            const core::DecompositionClass& decomposition_class) {
        CreateUniqueCDR3Map(decomposition_class);
        std::string cdrs_fasta = WriteUniqueCDR3InFasta(decomposition_class);
        std::string graph_fname = GetGraphFname(decomposition_class);
        return ComputeCDR3HammingGraphs(cdrs_fasta, graph_fname);
    }




}