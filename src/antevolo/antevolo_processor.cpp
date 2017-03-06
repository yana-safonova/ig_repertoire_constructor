#include <logger/logger.hpp>

#include "antevolo_processor.hpp"
#include "clone_set_decomposers/vj_clone_set_decomposer.hpp"
#include "vj_class_processors/vj_class_processor.hpp"


namespace antevolo {

    EvolutionaryTreeStorage AntEvoloProcessor::JoinEvolutionaryStoragesFromThreads() {
        EvolutionaryTreeStorage resulting_tree_storage ;
        for(auto it = thread_tree_storages_.begin(); it != thread_tree_storages_.end(); it++) {
            resulting_tree_storage.AppendArchive(*it);
        }
        auto const& tree = *resulting_tree_storage.cbegin();
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
        std::vector<size_t> fake_clone_indices(config_.run_params.num_threads);
        std::vector<size_t> reconstructed(config_.run_params.num_threads);
        std::vector<size_t> rejected(config_.run_params.num_threads);
        for (size_t i = 0; i < fake_clone_indices.size(); ++i) {
            fake_clone_indices[i] = (2 * i + 1) * total_number_of_reads_ ;
        }

#pragma omp parallel for schedule(dynamic)
        for(size_t i = 0; i < vj_decomposition.Size(); i++) {
            size_t thread_id = omp_get_thread_num();
            auto vj_class = vj_decomposition.GetClass(i);
            CloneSetWithFakesPtr fakes_clone_set_ptr(new CloneSetWithFakes(clone_set_));
            auto vj_class_processor = VJClassProcessor(fakes_clone_set_ptr,
                                                       config_.output_params,
                                                       config_.algorithm_params,
                                                       model,
                                                       clone_by_read_constructor_,
                                                       fake_clone_indices[thread_id],
                                                       reconstructed[thread_id],
                                                       rejected[thread_id]);
            vj_class_processor.CreateUniqueCDR3Map(vj_class);
            std::string cdrs_fasta = vj_class_processor.WriteUniqueCDR3InFasta(vj_class);
            std::string graph_fname = vj_class_processor.GetGraphFname(vj_class);
            TRACE("--------------------------");
            TRACE("CDR3 fasta: "<< cdrs_fasta << ", CDR3 Hamming graph: " << graph_fname);
            auto connected_components = vj_class_processor.ComputeCDR3HammingGraphs(cdrs_fasta, graph_fname);
            TRACE("# connected components: " << connected_components.size());
            for(size_t component_index = 0; component_index < connected_components.size(); component_index++) {
                EvolutionaryTree tree = vj_class_processor.AddComponent(
                        connected_components[component_index], component_index);
                tree.SetTreeIndices(i+1, component_index, 0);
                if (tree.NumEdges() != 0) {
                    thread_tree_storages_[thread_id].Add(tree);
                    //TRACE(i + 1 << "-th clonal tree was written to " << tree_output_fname);
                }
            }
        }
        size_t total_reconstructed = 0;
        size_t total_rejected = 0;
        for (size_t i = 0; i < reconstructed.size(); ++i) {
            total_reconstructed += reconstructed[i];
            total_rejected += rejected[i];
        }
        INFO("Number of reconstructed clones: " << total_reconstructed
             << ",\tnumber of edges rejected due to inequality of insertion blocks: " << total_rejected);
        return JoinEvolutionaryStoragesFromThreads();
    }
}