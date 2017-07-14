#pragma once

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "evolutionary_graph_utils/evolutionary_edge/base_evolutionary_edge.hpp"
#include "evolutionary_edge_constructor.hpp"
#include "../clone_set_with_fakes.hpp"

namespace antevolo {
    class EvolutionaryTree {
        CloneSetWithFakesPtr clone_set_ptr_;
        boost::unordered_map<size_t, EvolutionaryEdgePtr> edges_; // key is a src clone
        std::set<size_t> vertices_;

        std::vector<EvolutionaryEdgePtr> all_edge_vector_;
        std::map<size_t, std::vector<EvolutionaryEdgePtr>> outgoing_edges_;

        size_t VJ_class_index_;
        size_t connected_component_index_;
        size_t tree_index_;
        //std::string tree_output_fname_;
        //std::string vertices_output_fname_;

    public:
        EvolutionaryTree(CloneSetWithFakesPtr clone_set_ptr) : clone_set_ptr_(clone_set_ptr) {}

        void ReplaceEdge(size_t clone_num, EvolutionaryEdgePtr edge);

        void AddEdge(size_t dst_id, EvolutionaryEdgePtr edge);

        void AddDirected(size_t clone_num, EvolutionaryEdgePtr edge);

        void AddUndirected(size_t clone_num, EvolutionaryEdgePtr edge);

        void AddAllEdges();


        bool HasParentEdge(size_t clone_num) const;

        size_t NumEdges() const { return edges_.size(); }
        size_t NumVertices() const { return vertices_.size(); }

        typedef std::vector<EvolutionaryEdgePtr>::iterator EdgeIterator;

        typedef std::vector<EvolutionaryEdgePtr>::const_iterator ConstEdgeIterator;

        EdgeIterator begin() { return all_edge_vector_.begin(); }

        EdgeIterator end() { return all_edge_vector_.end(); }

        ConstEdgeIterator cbegin() const { return all_edge_vector_.cbegin(); }

        ConstEdgeIterator cend() const { return all_edge_vector_.cend(); }


        typedef std::set<size_t>::const_iterator ConstVertexIterator;

        ConstVertexIterator c_vertex_begin() const { return vertices_.cbegin(); }

        ConstVertexIterator c_vertex_end() const { return vertices_.cend(); }


        const EvolutionaryEdgePtr& GetParentEdge(size_t clone_num) const;

        // could be parent edge length or distance to the germline
        size_t GetParentEdgeLength(size_t clone_num) const;

        const std::vector<EvolutionaryEdgePtr>& OutgoingEdges(size_t clone_id) const;

        bool IsRoot(size_t clone_id) const;

        bool IsLeaf(size_t clone_id) const;

        bool IsFakeToFilter(size_t clone_id) const;

        size_t GetRoot() const;

        bool IsIsolated(size_t clone_id) const;

        bool ContainsClone(size_t clone_id) const {
            return vertices_.find(clone_id) != vertices_.end();
        }

        bool IsForest() const;

        size_t GetRootNumber() const;

        size_t GetRootByVertex(size_t clone_id) const;

        std::vector<size_t> GetRoots() const;

        size_t EdgeDepth() const;

        void SetTreeIndices(size_t VJ_class_index,
                            size_t connected_component_index,
                            size_t tree_index);

        std::string GetTreeOutputFname(std::string output_dir) const;

        size_t GetVJClassIndex() const { return VJ_class_index_;}
        size_t GetConnectedComponentIndex() const { return connected_component_index_;}
        size_t GetTreeIndex() const { return tree_index_;}

        const CloneSetWithFakes& GetCloneSet() const { return *clone_set_ptr_; }
        CloneSetWithFakesPtr GetCloneSetPtr() const { return clone_set_ptr_; }
    };

    std::ostream& operator<<(std::ostream& out, const EvolutionaryTree &tree);
}