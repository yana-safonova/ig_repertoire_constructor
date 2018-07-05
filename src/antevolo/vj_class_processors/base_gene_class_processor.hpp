#pragma once

#include "../clone_set_with_fakes.hpp"
#include <decomposition.hpp>
#include <antevolo_config.hpp>
#include "cdr3_hamming_graph_connected_components_processors/kruskal_cdr3_hg_cc_processor.hpp"
#include <shm_model_utils/shm_model_edge_weight_calculator.hpp>

namespace antevolo {
    class BaseGeneClassProcessor {
    protected:
        CloneSetWithFakesPtr clone_set_ptr_;

        const core::DecompositionClass &decomposition_class_;
        const AntEvoloConfig &config_;
        const AnnotatedCloneByReadConstructor &clone_by_read_constructor_;
        size_t current_fake_clone_index_;
        size_t reconstructed_;
        typedef std::map<std::string, std::vector<size_t>> UniqueCDR3IndexMap;
        typedef std::map<std::string, size_t> CDR3ToIndexMap;
        typedef boost::associative_property_map<std::map<size_t, size_t>> AP_map;

        UniqueCDR3IndexMap unique_cdr3s_map_;
        CDR3ToIndexMap cdr3_to_old_index_map_;
        std::vector<std::string> unique_cdr3s_;
        SparseGraphPtr sparse_cdr_graph_;
        GraphComponentMap graph_component_map_;


        void Clear();

        std::string GetFastaFname();

        virtual EvolutionaryTree ProcessComponentWithEdmonds(SparseGraphPtr hg_component, size_t component_id,
                                                             const ShmModelEdgeWeightCalculator &edge_weight_calculator);

        virtual void CreateUniqueCDR3Map() = 0;

        std::string WriteUniqueCDR3InFasta();

        std::string GetGraphFname();

        std::vector<SparseGraphPtr>
        ComputeCDR3HammingGraphs(std::string cdr_fasta, std::string cdr_graph, size_t tau);

    public:

        BaseGeneClassProcessor(CloneSetWithFakesPtr clone_set_ptr,
                               const core::DecompositionClass &decomposition_class,
                               const AntEvoloConfig &config,
                               const AnnotatedCloneByReadConstructor &clone_by_read_constructor,
                               size_t current_fake_clone_index) :
                clone_set_ptr_(clone_set_ptr),
                decomposition_class_(decomposition_class),
                config_(config),
                clone_by_read_constructor_(clone_by_read_constructor),
                current_fake_clone_index_(current_fake_clone_index),
                reconstructed_(0) {}

        virtual std::vector<SparseGraphPtr> ComputeConnectedComponents() = 0;

        virtual EvolutionaryTree ProcessComponent(SparseGraphPtr hg_component, size_t component_id,
                                                  const ShmModelEdgeWeightCalculator &edge_weight_calculator);

        size_t GetCurrentFakeCloneIndex() const { return current_fake_clone_index_; };

        size_t GetNumberOfReconstructedClones() const { return reconstructed_; };

        virtual ~BaseGeneClassProcessor() {}
    };

    typedef std::shared_ptr<BaseGeneClassProcessor> GeneCLassProcessorPtr;
}