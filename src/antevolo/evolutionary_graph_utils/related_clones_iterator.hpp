#pragma once

#include <cdr3_hamming_graph_connected_components_processors/base_cdr3_hg_cc_processor.hpp>
#include "../../graph_utils/sparse_graph.hpp"

namespace antevolo {

    enum RelatedClonesIteratorStage {SameSuperVertex, DifferentSuperVertices};

    using Base_CDR3_HG_CC_Processor::UniqueCDR3IndexMap;
    class RelatedClonesIterator {
        const SparseGraphPtr& hg_component_;
        size_t component_id_;
        size_t cdr3_index_;
        RelatedClonesIteratorStage stage_;
        size_t neighbour_super_vertex_num_;
        const UniqueCDR3IndexMap& unique_cdr3s_map_;
        const std::vector<std::string>& unique_cdr3s_;
        const GraphComponentMap& graph_component_map_;
        size_t related_clone_num_;

    public:
        RelatedClonesIterator(const SparseGraphPtr& hg_component,
                              size_t component_id,
                              size_t cdr3_index,
                              const UniqueCDR3IndexMap& unique_cdr3s_map,
                              const std::vector<std::string>& unique_cdr3s,
                              const GraphComponentMap& graph_component_map) :
                hg_component_(hg_component),
                component_id_(component_id),
                cdr3_index_(cdr3_index),
                unique_cdr3s_map_(unique_cdr3s_map),
                unique_cdr3s_(unique_cdr3s),
                graph_component_map_(graph_component_map),
                stage_(RelatedClonesIteratorStage::SameSuperVertex),
                neighbour_super_vertex_num_(0),
                related_clone_num_(0) {}

        RelatedClonesIterator& Next() {
            if (stage_ == RelatedClonesIteratorStage::SameSuperVertex) {
                const auto& clones_sharing_cdr3 = unique_cdr3s_map_.find(
                        unique_cdr3s_[neighbour_super_vertex_num_])->second;
                if (related_clone_num_ >= clones_sharing_cdr3.size()) {
                    stage_ = RelatedClonesIteratorStage::DifferentSuperVertices;
                    related_clone_num_ = 0;
                    neighbour_super_vertex_num_ = 0;
                }
            }
        }


    };

}
