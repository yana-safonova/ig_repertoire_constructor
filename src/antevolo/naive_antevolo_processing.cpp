#include <logger/logger.hpp>

#include "naive_antevolo_processing.hpp"
#include "clone_set_decomposers/vj_clone_set_decomposer.hpp"
#include "clonally_related_candidates_calculators/undirectred_first_tree_calculator.hpp"
#include "evolutionary_graph_utils/evolutionary_graph_constructor.hpp"


namespace antevolo {
    std::string NaiveAntEvoloProcessing::GetTreeOutputFname(std::string output_dir, size_t index1, size_t index2, size_t v_num, size_t e_num) {
        std::stringstream ss;
        ss << "clonal_tree_" << index1 << "-" << index2 << "_Vsize_" << v_num << "_Esize_" << e_num << ".tree";
        return path::append_path(output_dir, ss.str());
    }
    std::string NaiveAntEvoloProcessing::GetTreeClonesOutputFname(std::string output_dir, size_t index1, size_t index2, size_t v_num, size_t e_num) {
        std::stringstream ss;
        ss << "clonal_tree_" << index1 << "-" << index2 << "_Vsize_" << v_num << "_Esize_" << e_num << ".clones";
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
        ShmModel model(5, config_.input_params);
        for(size_t i = 0; i < vj_decomposition.Size(); i++) {
        //for(size_t i = 115; i < vj_decomposition.Size() && i < 116; i++) {
            auto vj_class = vj_decomposition.GetClass(i);
            auto candidate_calculator = UndirectedFirstTreeCalculator(clone_set_,
                                                                      config_.output_params,
                                                                      config_.algorithm_params,
                                                                      model);
            candidate_calculator.CreateUniqueCDR3Map(vj_class);
            std::string cdrs_fasta = candidate_calculator.WriteUniqueCDR3InFasta(vj_class);
            std::string graph_fname = candidate_calculator.GetGraphFname(vj_class);
            TRACE("--------------------------");
            TRACE("CDR3 fasta: "<< cdrs_fasta << ", CDR3 Hamming graph: " << graph_fname);
            auto connected_components = candidate_calculator.ComputeCDR3HammingGraphs(cdrs_fasta, graph_fname);
            TRACE("# connected components: " << connected_components.size());
            for(size_t component_index = 0; component_index < connected_components.size(); component_index++) {
            //for(size_t component_index = 2472; component_index < 2473 && component_index < connected_components.size(); component_index++) {
                EvolutionaryTree tree;
                candidate_calculator.AddComponent(
                        connected_components[component_index], component_index, tree);
                std::string tree_output_fname = GetTreeOutputFname(
                        config_.output_params.tree_dir, i + 1, component_index, tree.NumVertives(), tree.NumEdges());
                std::string vertices_output_fname = GetTreeOutputFname(
                        config_.output_params.vertex_dir, i + 1, component_index, tree.NumVertives(), tree.NumEdges());
                if (tree.NumEdges() != 0) {
                    tree.WriteInFileWithCDR3s(tree_output_fname);
                    tree.WriteVerticesInFile(vertices_output_fname, clone_set_);
                    TRACE(i + 1 << "-th clonal tree was written to " << tree_output_fname);
                }
            }
        }
        INFO("Clonal trees were written to " << config_.output_params.tree_dir);
    }
}