#include "base_cluster_to_forest_calculator.hpp"

namespace antevolo {
    class KruskalClusterToForestCalculator : public BaseClusterToForestCalculator {

        void SetUndirectedComponentParentEdge(size_t root_num, EvolutionaryEdgePtr edge);

        void PrepareSubtree(std::vector<std::pair<size_t, size_t>>& edge_vector,
                                                              size_t root_num);

        void PrepareSubtreeVertices(boost::unordered_set<size_t>& vertices_set,
                                    size_t root_num);

        void PrepareSubtreeKruskal(
                std::vector<std::pair<size_t, size_t>>& edge_vector,
                size_t root_vertex,
                const annotation_utils::CDRAnnotatedCloneSet& clone_set,
                std::shared_ptr<PolyEvolutionaryEdgeConstructor> edge_constructor);
    public:
        KruskalClusterToForestCalculator(const annotation_utils::CDRAnnotatedCloneSet &clone_set,
                                         const AntEvoloConfig::AlgorithmParams &config,
                                         GraphComponentMap& graph_component,
                                         const UniqueCDR3IndexMap& unique_cdr3s_map,
                                         const std::vector<std::string>& unique_cdr3s) :
                BaseClusterToForestCalculator(clone_set,
                                              config,
                                              graph_component,
                                              unique_cdr3s_map,
                                              unique_cdr3s) {}
        void SetUndirectedComponentsParentEdges(SparseGraphPtr hg_component,
                                                size_t component_id,
                                                boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges);
        void SetDirections(boost::unordered_set<size_t> vertices_nums,
                           EvolutionaryTree& tree,
                           boost::disjoint_sets<AP_map, AP_map> ds_on_undirected_edges);
    };
}