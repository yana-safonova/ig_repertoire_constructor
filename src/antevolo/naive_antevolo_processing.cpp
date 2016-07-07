#include <logger/logger.hpp>

#include "naive_antevolo_processing.hpp"
#include "clone_set_decomposers/vj_clone_set_decomposer.hpp"
#include "clonally_related_candidates_calculators/similar_cdr3_candidate_calculator.hpp"
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
            auto vj_class = vj_decomposition.GetClass(i);
            auto candidates_vector = SimilarCDR3CandidateCalculator(clone_set_,
                                                                    config_.output_params,
                                                                    config_.algorithm_params.similar_cdr3s_params.num_mismatches).
                    ComputeCandidates(vj_class);
            TRACE("Candidates vector contains " << candidates_vector.size() << " components");
            size_t candidate_index = 1;
            for(auto it = candidates_vector.begin(); it != candidates_vector.end(); it++) {
                if(it->Empty()) {
                    TRACE("Candidates are empty");
                    continue;
                }
                EvolutionaryGraphConstructor graph_constructor(config_.algorithm_params, clone_set_, *it);
                auto graph = graph_constructor.ConstructGraph();
                TRACE("Evolutionary graph size: " << graph.NumEdges());
                if(graph.NumEdges() == 0)
                    continue;
                auto tree = graph.GetEvolutionaryTree();
                std::string output_fname = GetTreeOutputFname(config_.output_params.tree_dir, i + 1, candidate_index,
                                                              tree.NumEdges());
                tree.WriteInFile(output_fname);
                TRACE(i + 1 << "-th clonal tree was written to " << output_fname);
                candidate_index++;
            }
        }
        INFO("Clonal trees were written to " << config_.output_params.tree_dir);
    }
}