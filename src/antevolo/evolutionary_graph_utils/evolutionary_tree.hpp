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
        std::set<size_t> vertices_;

        std::vector<EvolutionaryEdge> all_edge_vector_;
        std::map<size_t, std::vector<EvolutionaryEdge>> outgoing_edges_;

        size_t VJ_class_index_;
        size_t connected_component_index_;
        size_t tree_index_;
        //std::string tree_output_fname_;
        //std::string vertices_output_fname_;

        void AddEdge(size_t dst_id, EvolutionaryEdge edge);

    public:
        void AddDirected(size_t clone_num, EvolutionaryEdge edge);

        void AddUndirected(size_t clone_num, EvolutionaryEdge edge);

        bool Contains(size_t clone_num) const {
            return (edges_.find(clone_num) != edges_.end());
        }

        size_t NumEdges() const { return edges_.size(); }
        size_t NumVertices() const { return vertices_.size(); }

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

        void SetTreeIndices(size_t VJ_class_index,
                            size_t connected_component_index,
                            size_t tree_index);
        std::string GetTreeOutputFname(std::string output_dir) const;

        size_t GetVJClassIndex() const { return VJ_class_index_;}
        size_t GetConnectedComponentIndex() const { return connected_component_index_;}
        size_t GetTreeIndex() const { return tree_index_;}
    };

    std::ostream& operator<<(std::ostream& out, const EvolutionaryTree &tree);
}