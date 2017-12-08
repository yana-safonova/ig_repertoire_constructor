#pragma once

#include "base_cdr3_hg_cc_processor.hpp"
#include <boost/unordered_map.hpp>
#include <shm_model_utils/shm_model_edge_weight_calculator.hpp>
#include "edmonds_utils/edmonds_processor.hpp"

namespace antevolo {
    class Edmonds_CDR3_HG_CC_Processor : public Base_CDR3_HG_CC_Processor {
        boost::unordered_map<size_t, EvolutionaryEdgePtr> shorthest_directed_edge_;
        boost::unordered_set<size_t> vertices_nums_;
        const ShmModelEdgeWeightCalculator& edge_weight_calculator_;

        void SetShortestDirectedParentEdges();

        std::vector<WeightedEdge<int>> PrepareEdgeVector();

        void SetEdges(EvolutionaryTree& tree, const std::vector<WeightedEdge<int>>& edge_vector);


        double GetLength(EvolutionaryEdgePtr edge) {
            return edge_weight_calculator_.calculate_weigth_edge(*edge);
        }


    public:
        Edmonds_CDR3_HG_CC_Processor(CloneSetWithFakesPtr clone_set_ptr,
                                     const AntEvoloConfig::AlgorithmParams &config,
                                     const AnnotatedCloneByReadConstructor& clone_by_read_constructor,
                                     CDR3HammingGraphComponentInfo& hamming_graph_info,
                                     size_t current_fake_clone_index,
                                     const ShmModelEdgeWeightCalculator& edge_weight_calculator)
                : Base_CDR3_HG_CC_Processor(clone_set_ptr,
                                            config,
                                            clone_by_read_constructor,
                                            hamming_graph_info,
                                            current_fake_clone_index),
                   edge_weight_calculator_(edge_weight_calculator){}

        virtual EvolutionaryTree Process() override;
    };
}