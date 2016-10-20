#pragma once

#include "../../graph_utils/sparse_graph.hpp"
#include <evolutionary_graph_utils/evolutionary_tree.hpp>
#include <evolutionary_graph_utils/evolutionary_edge_constructor.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/unordered_set.hpp>

namespace antevolo {

    class BaseClusterToForestCalculator {

    public:

        typedef std::map<std::string, std::vector<size_t>> UniqueCDR3IndexMap;
        typedef boost::associative_property_map<std::map<size_t, size_t>> AP_map;

    protected:
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;
        const AntEvoloConfig::AlgorithmParams &config_;
        GraphComponentMap& graph_component_;
        const UniqueCDR3IndexMap& unique_cdr3s_map_;
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
        virtual void SetDirections(boost::unordered_set<size_t> vertices_nums,
                           EvolutionaryTree& tree,
                           boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges) = 0;
    public:

        BaseClusterToForestCalculator(const annotation_utils::CDRAnnotatedCloneSet &clone_set,
                                      const AntEvoloConfig::AlgorithmParams &config,
                                      GraphComponentMap& graph_component,
                                      const UniqueCDR3IndexMap& unique_cdr3s_map,
                                      const std::vector<std::string>& unique_cdr3s);

        EvolutionaryTree ConstructForest(SparseGraphPtr hg_component, size_t component_id);

        std::shared_ptr<PolyEvolutionaryEdgeConstructor> GetEdgeConstructor() {
            PolyEvolutionaryEdgeConstructor* ptr = new PolyVJEvolutionaryEdgeConstructor(config_.edge_construction_params);
            return std::shared_ptr<PolyEvolutionaryEdgeConstructor>(ptr);
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

        virtual ~BaseClusterToForestCalculator() {};
    };


}