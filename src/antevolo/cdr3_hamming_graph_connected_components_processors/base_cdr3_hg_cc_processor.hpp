#pragma once

#include "../../graph_utils/sparse_graph.hpp"
#include <evolutionary_graph_utils/evolutionary_edge_constructor.hpp>
#include <evolutionary_graph_utils/evolutionary_tree.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/unordered_set.hpp>
#include <annotated_clone_by_read_constructor.hpp>
#include "../clone_set_with_fakes.hpp"
#include <evolutionary_graph_utils/related_clones_iterator.hpp>

namespace antevolo {

    class Base_CDR3_HG_CC_Processor {

    public:
        typedef std::map<std::string, std::vector<size_t>> UniqueCDR3IndexMap;
        typedef std::map<std::string, size_t> CDR3ToIndexMap;
        typedef boost::associative_property_map<std::map<size_t, size_t>> AP_map;

    protected:
//        const annotation_utils::CDRAnnotatedCloneSet& clone_set_;
        CloneSetWithFakesPtr clone_set_ptr_;
        const AntEvoloConfig::AlgorithmParams& config_;
        const AnnotatedCloneByReadConstructor& clone_by_read_constructor_;
        CDR3HammingGraphComponentInfo& hamming_graph_info_;
        size_t current_fake_clone_index_;
        size_t reconstructed_;

        boost::unordered_map<size_t, std::set<size_t>> undirected_graph_;
        boost::unordered_map<size_t, bool> parent_edge_handled_;
        boost::unordered_map<size_t, EvolutionaryEdgePtr> undirected_components_edges_;

        static const size_t EVO_EDGE_MAX_LENGTH = 400; // todo: move to config

        void AddUndirectedPair(size_t src_num, size_t dst_num);

        void AddUndirectedForest(boost::disjoint_sets<AP_map, AP_map> &ds_on_undirected_edges,
                                 const boost::unordered_set<size_t>& vertices_nums);

//        virtual void SetUndirectedComponentsParentEdges(boost::disjoint_sets<AP_map, AP_map>& ds_on_undirected_edges,
//                                                        const boost::unordered_set<size_t>& vertices_nums) = 0;
//
//        virtual void SetDirections(boost::disjoint_sets<AP_map, AP_map>& ds_on_undirected_edges,
//                                   const boost::unordered_set<size_t> &vertices_nums,
//                                   EvolutionaryTree &tree) = 0;

        void ReconstructMissingVertices(boost::unordered_set<size_t>& vertices_nums,
                                                EvolutionaryTree& tree);

        void Refine(boost::unordered_set<size_t>& vertices_nums,
                    EvolutionaryTree& tree);

        bool SecondCloneIsFirstsAncestor(EvolutionaryTree& tree, size_t first_clone, size_t second_clone);


        bool ReconstructAncestralLineageSimple(
                EvolutionaryEdgePtr edge,
                EvolutionaryTree &tree,
                boost::unordered_set<size_t> &vertices_nums,
                const std::shared_ptr<EvolutionaryEdgeConstructor> &edge_constructor,
                std::vector<size_t> &roots,
                boost::unordered_map<size_t, size_t> &iterator_index_map);
        void HandleRootNeighbour(
                size_t root_num,
                size_t dst_num,
                boost::unordered_set<size_t>& vertices_nums,
                EvolutionaryTree& tree,
                boost::unordered_map<size_t, EvolutionaryEdgePtr>& roots_nearest_neighbours,
                const std::shared_ptr<EvolutionaryEdgeConstructor>& edge_constructor);

        size_t GetUndirectedCompopentRoot(size_t root_num) {
            if (undirected_components_edges_.find(root_num) != undirected_components_edges_.end()) {
                return undirected_components_edges_[root_num]->DstNum();
            }
            return size_t(-1);
        }

        const EvolutionaryEdgePtr& GetUndirectedComponentParentEdge(size_t root_num) {
            return undirected_components_edges_[root_num];
        }

    public:

        Base_CDR3_HG_CC_Processor(CloneSetWithFakesPtr clone_set_ptr,
                                  const AntEvoloConfig::AlgorithmParams &config,
                                  const AnnotatedCloneByReadConstructor& clone_by_read_constructor,
                                  CDR3HammingGraphComponentInfo& hamming_graph_info,
                                  size_t current_fake_clone_index);

        virtual EvolutionaryTree Process() = 0;

        std::shared_ptr<EvolutionaryEdgeConstructor> GetEdgeConstructor() {
            EvolutionaryEdgeConstructor* ptr = new VJEvolutionaryEdgeConstructor(config_.edge_construction_params);
            return std::shared_ptr<EvolutionaryEdgeConstructor>(ptr);
        }

        static bool CheckClonesConsistencyForReconstruction(const annotation_utils::AnnotatedClone& left,
                                                            const annotation_utils::AnnotatedClone& right);

        size_t GetCurrentFakeCloneIndex() const { return current_fake_clone_index_; };

        size_t GetNumberOfReconstructedClones() const { return reconstructed_; };

        virtual ~Base_CDR3_HG_CC_Processor() {};
    };


}