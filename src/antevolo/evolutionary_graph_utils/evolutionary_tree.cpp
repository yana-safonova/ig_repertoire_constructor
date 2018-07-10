#include <verify.hpp>
#include "evolutionary_tree.hpp"

namespace antevolo {
    void EvolutionaryTree::AddEdge(size_t dst_id, EvolutionaryEdgePtr edge) {
        VERIFY(dst_id == edge->DstNum());
        edges_[dst_id] = edge;
        all_edge_vector_.push_back(edge);
        vertices_.insert(edge->DstNum());
        vertices_.insert(edge->SrcNum());
        if(outgoing_edges_.find(edge->SrcNum()) == outgoing_edges_.end()) {
            outgoing_edges_[edge->SrcNum()] = std::vector<EvolutionaryEdgePtr>();
        }
        outgoing_edges_[edge->SrcNum()].push_back(edge);
    }

    void EvolutionaryTree::AddDirected(size_t clone_num, EvolutionaryEdgePtr edge) {
        VERIFY(edge->IsDirected());
        //std::cout << "oppa: " << clone_num << " - " << edge.src_clone_num << " - " << edge.dst_clone_num << std::endl;
        if(edge->IsDirected()) {
            if (!HasParentEdge(clone_num)) {
                ReplaceEdge(clone_num, edge);
                return;
            }
            const EvolutionaryEdgePtr& parent_edge =  edges_[clone_num];
            if (parent_edge->Length() > edge->Length()) { //todo: compare only num added shms ?
                //if clone_set_[*it2] is root or if the new edge is shorter
                ReplaceEdge(clone_num, edge);
                return;
            }
            /*
            if (parent_edge.num_added_shms == edge.num_added_shms && parent_edge.cdr3_distance > edge.cdr3_distance) {
                //if clone_set_[*it2] is root or if the new edge is shorter
                AddEdge(clone_num, edge);
                return;
            }
            */
        }
    }

    void EvolutionaryTree::AddUndirected(size_t clone_num, EvolutionaryEdgePtr edge) {
        VERIFY(edge->IsUndirected());
        if (edge->IsUndirected()) {
            ReplaceEdge(clone_num, edge);
        }
    }

    void EvolutionaryTree::ReplaceEdge(size_t clone_num, EvolutionaryEdgePtr edge) {
        VERIFY(edge->DstNum() == clone_num);
        edges_[clone_num] = edge;
    }

    void EvolutionaryTree::AddAllEdges() {
        for (auto p : edges_) {
            VERIFY(p.second->DstClone()->CDR3Range().length() == p.second->SrcClone()->CDR3Range().length());
            AddEdge(p.first, p.second);
        }
    }

    bool EvolutionaryTree::HasParentEdge(size_t clone_num) const {
        return edges_.find(clone_num) != edges_.end();
    }

    void EvolutionaryTree::AddVertex(size_t clone_num) {
        vertices_.insert(clone_num);
    }

    const EvolutionaryEdgePtr& EvolutionaryTree::GetParentEdge(size_t clone_num) const {
        auto p = edges_.find(clone_num);
        VERIFY_MSG(p != edges_.end(), "evolutionary tree: got a request for unexisting edge");
        return p->second;
    }

    size_t EvolutionaryTree::GetParentEdgeLength(size_t clone_num) const {
        if (HasParentEdge(clone_num)) {
            return GetParentEdge(clone_num)->Length();
        }
        else {
            const auto& clone = clone_set_ptr_->operator[](clone_num);
            return clone.VSHMs().size() + clone.JSHMs().size();
        }
    }

    bool EvolutionaryTree::IsRoot(size_t clone_id) const {
        VERIFY_MSG(vertices_.find(clone_id) != vertices_.end(), "Tree does not contain vertex " << clone_id);
        return edges_.find(clone_id) == edges_.end();
    }

    size_t EvolutionaryTree::GetRoot() const {
        for(auto it = vertices_.begin(); it != vertices_.end(); it++) {
            if(IsRoot(*it)) {
                return *it;
            }
        }
        VERIFY_MSG(false, "Root was not found");
        return size_t(-1);
    }

    bool EvolutionaryTree::IsLeaf(size_t clone_id) const {
        VERIFY_MSG(ContainsClone(clone_id), "Tree does not contain vertex " << clone_id);
        return outgoing_edges_.find(clone_id) == outgoing_edges_.end();
    }

    bool EvolutionaryTree::IsFakeToFilter(size_t clone_id) const {
        VERIFY_MSG(ContainsClone(clone_id), "Tree does not contain vertex " << clone_id);
        return clone_set_ptr_->IsFake(clone_id) && (IsLeaf(clone_id) || OutgoingEdges(clone_id).size() < 2);
    }

    bool EvolutionaryTree::IsForest() const {
        size_t num_roots = 0;
        for(auto it = vertices_.begin(); it != vertices_.end(); it++) {
            if(IsRoot(*it)) {
                num_roots++;
            }
        }
        return num_roots > 1;
    }

    size_t EvolutionaryTree::GetRootNumber() const {
        size_t num_roots = 0;
        for(auto it = vertices_.begin(); it != vertices_.end(); it++) {
            if(IsRoot(*it)) {
                num_roots++;
            }
        }
        return num_roots;
    }

    size_t EvolutionaryTree::EdgeDepth() const {
        VERIFY_MSG(false, "Implement me!");
        return 0;
    }

    const std::vector<EvolutionaryEdgePtr>& EvolutionaryTree::OutgoingEdges(size_t clone_id) const {
        VERIFY_MSG(ContainsClone(clone_id), "Tree does not contain vertex " << clone_id);
        VERIFY_MSG(!IsLeaf(clone_id), "Vertex " << clone_id << " is leaf");
        const std::vector<EvolutionaryEdgePtr>& v = outgoing_edges_.at(clone_id);
        return v; // hmm?
    }

    std::vector<size_t> EvolutionaryTree::GetRoots() const {
        std::vector<size_t> roots;
        for(auto it = vertices_.begin(); it != vertices_.end(); it++) {
            if(IsRoot(*it)) {
                roots.push_back(*it);
            }
        }
        return roots;
    }

    void EvolutionaryTree::SetTreeIndices(size_t VJ_class_index,
                                          size_t connected_component_index,
                                          size_t tree_index) {
        VJ_class_index_ = VJ_class_index;
        connected_component_index_ = connected_component_index;
        tree_index_ = tree_index;
    }
    std::string EvolutionaryTree::GetTreeOutputFname(std::string output_dir) const{
        std::stringstream ss;
        ss << "clonal_tree_" << VJ_class_index_ << "-" << connected_component_index_ << "-" << tree_index_
        << "_Vsize_" << NumVertices() << "_Esize_" << NumEdges() << ".tree";
        return path::append_path(output_dir, ss.str());
    }

    bool EvolutionaryTree::IsIsolated(size_t clone_id) const {
        return IsRoot(clone_id) and IsLeaf(clone_id);
    }

    size_t EvolutionaryTree::GetRootByVertex(size_t clone_id) const {
        if(IsRoot(clone_id))
            return clone_id;
        size_t cur_parent = clone_id;
        do {
            cur_parent = GetParentEdge(cur_parent)->SrcNum();
        }
        while(!IsRoot(cur_parent));
        return cur_parent;
    }

    std::ostream& operator<<(std::ostream& out, const EvolutionaryTree &tree) {
        for(auto it = tree.cbegin(); it != tree.cend(); it++)
            out << *it << std::endl;
        return out;
    }
}
