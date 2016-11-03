#pragma once

#include "../../graph_utils/sparse_graph.hpp"
#include <evolutionary_graph_utils/evolutionary_edge_constructor.hpp>
#include <evolutionary_graph_utils/evolutionary_tree.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/unordered_set.hpp>
#include "../clone_set_with_fakes.hpp"

namespace antevolo {

    class Base_CDR3_HG_CC_Processor {

    public:

        typedef std::map<std::string, std::vector<size_t>> UniqueCDR3IndexMap;
        typedef std::map<std::string, size_t> CDR3ToIndexMap;
        typedef boost::associative_property_map<std::map<size_t, size_t>> AP_map;

    protected:
        CloneSetWithFakes& clone_set_;
        const AntEvoloConfig::AlgorithmParams &config_;
        GraphComponentMap& graph_component_;
        const UniqueCDR3IndexMap& unique_cdr3s_map_;
        const CDR3ToIndexMap& cdr3_to_index_map_;
        const std::vector<std::string>& unique_cdr3s_;

        boost::unordered_map<size_t, std::set<size_t>> undirected_graph_;
        boost::unordered_map<size_t, bool> parent_edge_handled_;
        boost::unordered_map<size_t, EvolutionaryEdgePtr> undirected_components_edges_;


        void AddUndirectedPair(size_t src_num, size_t dst_num);

        void AddUndirectedForest(SparseGraphPtr hg_component, size_t component_id,
                                 boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges);
        virtual void SetUndirectedComponentsParentEdges(SparseGraphPtr hg_component,
                                                size_t component_id,
                                                boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges) = 0;
        virtual void SetDirections(const boost::unordered_set<size_t>& vertices_nums,
                                   EvolutionaryTree& tree,
                                   boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges) = 0;
        virtual void ReconstructMissingVertices(const boost::unordered_set<size_t> &vertices_nums,
                                                        EvolutionaryTree &tree, SparseGraphPtr hg_component,
                                                        size_t component_id) = 0;
    public:

        Base_CDR3_HG_CC_Processor(CloneSetWithFakes& clone_set,
                                  const AntEvoloConfig::AlgorithmParams &config,
                                  GraphComponentMap& graph_component,
                                  const UniqueCDR3IndexMap& unique_cdr3s_map,
                                  const CDR3ToIndexMap& cdr3_to_index_map,
                                  const std::vector<std::string>& unique_cdr3s);

        EvolutionaryTree ConstructForest(SparseGraphPtr hg_component, size_t component_id);

        std::shared_ptr<EvolutionaryEdgeConstructor> GetEdgeConstructor() {
            EvolutionaryEdgeConstructor* ptr = new VJEvolutionaryEdgeConstructor(config_.edge_construction_params);
            return std::shared_ptr<EvolutionaryEdgeConstructor>(ptr);
        }

        size_t GetUndirectedCompopentRoot(size_t root_num) {
            if (undirected_components_edges_.find(root_num) != undirected_components_edges_.end()) {
                return undirected_components_edges_[root_num]->DstNum();
            }
            return size_t(-1);
        }

        const EvolutionaryEdgePtr& GetUndirectedComponentParentEdge(size_t root_num) {
            return undirected_components_edges_[root_num];
        }

        virtual ~Base_CDR3_HG_CC_Processor() {};
    };


}