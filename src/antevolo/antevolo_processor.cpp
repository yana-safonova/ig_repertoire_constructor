#include <logger/logger.hpp>

#include "antevolo_processor.hpp"
#include "clone_set_decomposers/vj_clone_set_decomposer.hpp"
#include "clonally_related_candidates_calculators/undirectred_first_tree_calculator.hpp"
#include "evolutionary_graph_utils/evolutionary_graph_constructor.hpp"


namespace antevolo {
    std::string AntEvoloProcessor::GetTreeOutputFname(std::string output_dir, size_t index1, size_t index2, size_t v_num, size_t e_num) {
        std::stringstream ss;
        ss << "clonal_tree_" << index1 << "-" << index2 << "_Vsize_" << v_num << "_Esize_" << e_num << ".tree";
        return path::append_path(output_dir, ss.str());
    }
    std::string AntEvoloProcessor::GetTreeClonesOutputFname(std::string output_dir, size_t index1, size_t index2, size_t v_num, size_t e_num) {
        std::stringstream ss;
        ss << "clonal_tree_" << index1 << "-" << index2 << "_Vsize_" << v_num << "_Esize_" << e_num << ".clones";
        return path::append_path(output_dir, ss.str());
    }

    EvolutionaryTreeStorage AntEvoloProcessor::JoinEvolutionaryStoragesFromThreads() {
        EvolutionaryTreeStorage resulting_tree_storage(clone_set_)  ;
        for(auto it = thread_tree_storages_.begin(); it != thread_tree_storages_.end(); it++)
            resulting_tree_storage.AppendArchive(*it);
        return resulting_tree_storage;
    }

    EvolutionaryTreeStorage AntEvoloProcessor::ConstructClonalTrees() {
        VJCloneSetDecomposer clone_set_decomposer(clone_set_);
        auto vj_decomposition = clone_set_decomposer.CreateDecomposition();
        INFO("VJ decomposition containing " << vj_decomposition.Size() << " classes was created.");
        INFO("Largest class contains " << vj_decomposition.MaxClassSize() << " clone(s)");
        omp_set_num_threads(config_.run_params.num_threads);
        INFO("Construction of clonal trees starts");
        ShmModel model(5, config_.input_params); // todo: move magic constant to config
#pragma omp parallel for schedule(dynamic)
        for(size_t i = 0; i < vj_decomposition.Size(); i++) {
            size_t thread_id = omp_get_thread_num();
        //for(size_t i = 59; i < vj_decomposition.Size() && i < 60; i++) {
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
            //for(size_t component_index = 1; component_index < 2 && component_index < connected_components.size(); component_index++) {
                // todo: do not pass resulting evolutionary tree as a method argument, use return
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
                    thread_tree_storages_[thread_id].Add(tree);
                    TRACE(i + 1 << "-th clonal tree was written to " << tree_output_fname);
                }
            }
        }
        INFO("Clonal trees were written to " << config_.output_params.tree_dir);
        return JoinEvolutionaryStoragesFromThreads();
    }
}