#include <logger/logger.hpp>

#include "antevolo_processor.hpp"
#include "clone_set_decomposers/vj_clone_set_decomposer.hpp"
#include "vj_class_processors/vj_class_processor.hpp"


namespace antevolo {

    EvolutionaryTreeStorage AntEvoloProcessor::JoinEvolutionaryStoragesFromThreads() {
        EvolutionaryTreeStorage resulting_tree_storage(clone_set_)  ;
        for(auto it = thread_tree_storages_.begin(); it != thread_tree_storages_.end(); it++)
            resulting_tree_storage.AppendArchive(*it);
        return resulting_tree_storage;
    }

    EvolutionaryTreeStorage AntEvoloProcessor::ConstructClonalTrees() {

        VJCloneSetDecomposer clone_set_decomposer(clone_set_); // storage for reconstructed fake vertices
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
            auto candidate_calculator = VJClassProcessor(clone_set_,
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
                EvolutionaryTree tree = candidate_calculator.AddComponent(
                        connected_components[component_index], component_index);
                tree.SetTreeIndices(i+1, component_index, 0);
                if (tree.NumEdges() != 0) {
                    thread_tree_storages_[thread_id].Add(tree);
                    //TRACE(i + 1 << "-th clonal tree was written to " << tree_output_fname);
                }
            }
        }
        return JoinEvolutionaryStoragesFromThreads();
    }
}