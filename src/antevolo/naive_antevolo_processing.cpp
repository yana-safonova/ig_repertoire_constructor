#include <logger/logger.hpp>

#include "naive_antevolo_processing.hpp"
#include "clone_set_decomposers/vj_clone_set_decomposer.hpp"
#include "clonally_related_candidates_calculators/undirectred_first_tree_calculator.hpp"
#include "evolutionary_graph_utils/evolutionary_graph_constructor.hpp"

namespace antevolo {
    std::string NaiveAntEvoloProcessing::GetTreeOutputFname(std::string output_dir, size_t index1, size_t index2, size_t tree_size) {
        std::stringstream ss;
        ss << "clonal_tree_" << index1 << "-" << index2 << "_size_" << tree_size << ".tree";
        return path::append_path(output_dir, ss.str());
    }

    void NaiveAntEvoloProcessing::ConstructClonalTrees() {
        VJCloneSetDecomposer clone_set_decomposer(clone_set_);
        auto vj_decomposition = clone_set_decomposer.CreateDecomposition();
        INFO("VJ decomposition containing " << vj_decomposition.Size() << " classes was created.");
        INFO("Largest class contains " << vj_decomposition.MaxClassSize() << " clone(s)");
        //omp_set_num_threads(16); // todo: add in config
        INFO("Construction of clonal trees starts");
//#pragma omp parallel for schedule(dynamic)
        for(size_t i = 0; i < vj_decomposition.Size(); i++) {
        //for(size_t i = 301; i < 302 && i < vj_decomposition.Size(); i++) {
            auto vj_class = vj_decomposition.GetClass(i);
            auto candidate_calculator = UndirectedFirstTreeCalculator(clone_set_,
                                                                    config_.output_params,
                                                                    config_.algorithm_params);
            candidate_calculator.CreateUniqueCDR3Map(vj_class);
            std::string cdrs_fasta = candidate_calculator.WriteUniqueCDR3InFasta(vj_class);
            std::string graph_fname = candidate_calculator.GetGraphFname(vj_class);
            TRACE("--------------------------");
            TRACE("CDR3 fasta: "<< cdrs_fasta << ", CDR3 Hamming graph: " << graph_fname);
            auto connected_components = candidate_calculator.ComputeCDR3HammingGraphs(cdrs_fasta, graph_fname);
            TRACE("# connected components: " << connected_components.size());
            for(size_t component_index = 0; component_index < connected_components.size(); component_index++) {
            //for(size_t component_index = 5; component_index < 6 && component_index < connected_components.size(); component_index++) {
                EvolutionaryTree tree;
                candidate_calculator.AddUndirectedForestToTheTree(
                        connected_components[component_index], component_index, tree);
                candidate_calculator.AddComponentToTheTree(
                        connected_components[component_index], component_index, tree);
                std::string output_fname = GetTreeOutputFname(
                        config_.output_params.tree_dir, i + 1, component_index, tree.NumEdges());
                if (tree.NumEdges() != 0) {
                    tree.WriteInFile(output_fname);
                    TRACE(i + 1 << "-th clonal tree was written to " << output_fname);
                }
            }
        }
        INFO("Clonal trees were written to " << config_.output_params.tree_dir);
    }
}