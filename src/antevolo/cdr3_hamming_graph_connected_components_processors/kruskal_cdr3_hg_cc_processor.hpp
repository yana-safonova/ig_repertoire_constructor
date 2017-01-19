#include "base_cdr3_hg_cc_processor.hpp"

namespace antevolo {
    class Kruskal_CDR3_HG_CC_Processor : public Base_CDR3_HG_CC_Processor {

        void SetUndirectedComponentParentEdge(size_t root_num, EvolutionaryEdgePtr edge);

        void PrepareSubtree(std::vector<std::pair<size_t, size_t>>& edge_vector,
                                                              size_t root_num);

        void PrepareSubtreeVertices(boost::unordered_set<size_t>& vertices_set,
                                    size_t root_num);

        void PrepareSubtreeKruskal(
                std::vector<std::pair<size_t, size_t>>& edge_vector,
                size_t root_vertex,
                CloneSetWithFakes& clone_set,
                std::shared_ptr<EvolutionaryEdgeConstructor> edge_constructor);
        void ReconstructAncestralLineageSimple(
                EvolutionaryEdgePtr edge,
                EvolutionaryTree& tree,
                const std::shared_ptr<EvolutionaryEdgeConstructor>& edge_constructor,
                std::vector<size_t>& roots,
                boost::unordered_map<size_t, size_t>& iterator_index_map);
//        void ReconstructAncestralLineage(EvolutionaryEdgePtr edge,
//                                         EvolutionaryTree& tree,
//                                         const std::shared_ptr<EvolutionaryEdgeConstructor> &edge_constructor,
//                                         std::vector<size_t>& roots);

        void SetUndirectedComponentsParentEdges(boost::disjoint_sets<AP_map, AP_map>& ds_on_undirected_edges,
                                                const boost::unordered_set<size_t>& vertices_nums) override;
        void SetDirections(boost::disjoint_sets<AP_map, AP_map>& ds_on_undirected_edges,
                                   const boost::unordered_set<size_t> &vertices_nums, EvolutionaryTree &tree) override;

        void ReconstructMissingVertices(const boost::unordered_set<size_t> &vertices_nums,
                                                EvolutionaryTree &tree) override;

        void HandleRootNeighbour(
                size_t root_num,
                size_t dst_num,
                const boost::unordered_set<size_t>& vertices_nums,
                EvolutionaryTree& tree,
                boost::unordered_map<size_t, EvolutionaryEdgePtr>& roots_nearest_neighbours,
                const std::shared_ptr<EvolutionaryEdgeConstructor>& edge_constructor);
    public:
        Kruskal_CDR3_HG_CC_Processor(CloneSetWithFakes& clone_set,
                                     const AntEvoloConfig::AlgorithmParams &config,
                                     const AnnotatedCloneByReadConstructor& clone_by_read_constructor,
                                     CDR3HammingGraphInfo& hamming_graph_info)
                : Base_CDR3_HG_CC_Processor(clone_set,
                                            config,
                                            clone_by_read_constructor,
                                            hamming_graph_info) {}


    };
}