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
                std::shared_ptr<EvolutionaryEdgeConstructor> edge_constructor);

        void SetUndirectedComponentsParentEdges(boost::disjoint_sets<AP_map, AP_map>& ds_on_undirected_edges,
                                                const boost::unordered_set<size_t>& vertices_nums);

        void SetDirections(boost::disjoint_sets<AP_map, AP_map>& ds_on_undirected_edges,
                                   const boost::unordered_set<size_t> &vertices_nums, EvolutionaryTree &tree);
    public:

        EvolutionaryTree Process() override;

        Kruskal_CDR3_HG_CC_Processor(CloneSetWithFakesPtr clone_set_ptr,
                                     const AntEvoloConfig::AlgorithmParams &config,
                                     const AnnotatedCloneByReadConstructor& clone_by_read_constructor,
                                     CDR3HammingGraphComponentInfo& hamming_graph_info,
                                     size_t current_fake_clone_index)
                : Base_CDR3_HG_CC_Processor(clone_set_ptr,
                                            config,
                                            clone_by_read_constructor,
                                            hamming_graph_info,
                                            current_fake_clone_index) {}


    };
}