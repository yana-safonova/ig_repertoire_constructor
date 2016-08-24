#pragma once

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "model_utils/shm_model.hpp"
#include "evolutionary_edge.hpp"
#include "evolutionary_edge_constructor.hpp"
#include <annotation_utils/annotated_clone_set.hpp>

namespace antevolo {
    class EvolutionaryTree {
        boost::unordered_map<size_t, EvolutionaryEdge> edges_; // key is a src clone
        boost::unordered_map<size_t, std::set<size_t>> undirected_graph_;
        boost::unordered_map<size_t, bool> passed_flag_;
        boost::unordered_map<size_t, EvolutionaryEdge> undirected_components_edges_;
        std::set<size_t> vertices_;

        std::vector<EvolutionaryEdge> all_edge_vector_;
        std::map<size_t, std::vector<EvolutionaryEdge>> outgoing_edges_;

        void AddEdge(size_t dst_id, EvolutionaryEdge edge);

    public:
        void AddDirected(size_t clone_num, EvolutionaryEdge edge);

        void SetUndirectedComponentParentEdge(size_t root_num, EvolutionaryEdge edge, ShmModel& model);

        void AddUndirected(size_t clone_num, EvolutionaryEdge edge);

        void AddUndirectedPair(size_t src_num, size_t dst_num);

        void PrepareSubtree(std::vector<std::pair<size_t, size_t>>& edge_vector, size_t root_num);

        void PrepareSubtreeVertices(
                boost::unordered_set<size_t>& vertices_set,
                size_t root_num);

        void PrepareSubtreeKruskal(std::vector<std::pair<size_t, size_t>>& edge_vector,
                                                     size_t root_vertex,
                                                     const annotation_utils::CDRAnnotatedCloneSet& clone_set,
                                                     EvolutionaryEdgeConstructor* edge_constructor);

        void PrepareSubtreeEdmonds(std::vector<std::pair<size_t, size_t>>& edge_vector,
                                                     size_t root_vertex,
                                                     ShmModel& model,
                                                     const annotation_utils::CDRAnnotatedCloneSet& clone_set,
                                                     EvolutionaryEdgeConstructor* edge_constructor);

        boost::unordered_map<size_t, std::set<size_t>>& GetUndirectedGraph() {
            return undirected_graph_;
        };

        bool Contains(size_t clone_num) {
            return (edges_.find(clone_num) != edges_.end());
        }

        size_t GetUndirectedCompopentRoot(size_t root_num) {
            if (undirected_components_edges_.find(root_num) != edges_.end()){
                return undirected_components_edges_[root_num].dst_clone_num;
            }
            else {
                return size_t(-1);
            }
        }

        const EvolutionaryEdge& GetUndirectedComponentParentEdge(size_t root_num) {
            return undirected_components_edges_[root_num];
        }

        // todo: move all output methods from EvolutionaryTree
        void WriteEdge(const EvolutionaryEdge& edge, std::ofstream& out); //no endl

        void WriteInFile(std::string output_fname);
        void WriteInFileWithCDR3s(std::string output_fname);

        void WriteVerticesInFile(std::string output_fname, const annotation_utils::CDRAnnotatedCloneSet& clone_set);

        size_t NumEdges() const { return edges_.size(); }
        size_t NumVertices() const { return passed_flag_.size(); }

        size_t UndirectedGraphSize() const {
            size_t res = 0;
            for (auto it = undirected_graph_.begin(); it != undirected_graph_.end(); it++) {
                res += it->second.size();
            }
            return res;
        }

        // what are these methods about?
        void SetFlag(bool b, size_t clone_num) {
            passed_flag_[clone_num] = b;
        }
        bool GetFlag(size_t clone_num) {
            return passed_flag_[clone_num];
        }

        // why clone to string is a part of EvolutionaryTree class?
        std::string clone_to_string(const annotation_utils::AnnotatedClone& clone) {
            std::stringstream ss;
            size_t start_pos = clone.GetRangeByRegion(
                    annotation_utils::StructuralRegion::CDR3).start_pos;
            size_t end_pos = clone.GetRangeByRegion(
                    annotation_utils::StructuralRegion::CDR3).end_pos;
            auto const& seq = clone.Read().seq;
            for (size_t i = 0; i < seqan::length(seq); ++i) {
                if (i == start_pos) {
                    ss << " ";
                }
                ss << seq[i];
                if (i == end_pos) {
                    ss << " ";
                }
            }
            return ss.str();
        }

        typedef std::vector<EvolutionaryEdge>::iterator EdgeIterator;

        typedef std::vector<EvolutionaryEdge>::const_iterator ConstEdgeIterator;

        EdgeIterator begin() { return all_edge_vector_.begin(); }

        EdgeIterator end() { return all_edge_vector_.end(); }

        ConstEdgeIterator cbegin() const { return all_edge_vector_.cbegin(); }

        ConstEdgeIterator cend() const { return all_edge_vector_.cend(); }


        typedef std::set<size_t>::const_iterator ConstVertexIterator;

        ConstVertexIterator c_vertex_begin() const { return vertices_.cbegin(); }

        ConstVertexIterator c_vertex_end() const { return vertices_.cend(); }


        const EvolutionaryEdge& GetParentEdge(size_t clone_num) const {
            return edges_.at(clone_num);
        }

        const std::vector<EvolutionaryEdge>& OutgoingEdges(size_t clone_id) const;


        bool IsRoot(size_t clone_id) const;

        bool IsLeaf(size_t clone_id) const;

        size_t GetRoot() const;

        bool ContainsClone(size_t clone_id) const {
            return vertices_.find(clone_id) != vertices_.end();
        }

        bool IsForest() const;

        size_t GetRootNumber() const;

        std::vector<size_t> GetRoots() const;

        size_t EdgeDepth() const;
    };

    std::ostream& operator<<(std::ostream& out, const EvolutionaryTree &tree);
}