#include "clonal_graph_constructor.hpp"

namespace antevolo {
    void ClonalGraphConstructor::AddEdgesFromTree() {
        for(auto it = tree_.cbegin(); it != tree_.cend(); it++) {
            clonal_graph_.AddOldEdge((*it)->SrcNum(), (*it)->DstNum());
        }
    }

    bool ClonalGraphConstructor::SkipVertex(size_t vertex) const {
        if(tree_.IsRoot(vertex))
            return false;
        return tree_.GetParentEdge(vertex)->IsUndirected();
    }

    bool ClonalGraphConstructor::VerticesConnectedByDirectedPath(size_t v_src, size_t v_dst) const {
        if(tree_.IsRoot(v_dst))
            return false;
        size_t cur_parent = v_dst;
        do {
            auto parent_edge = tree_.GetParentEdge(cur_parent);
            //if(!parent_edge->IsDirected())
            //    return false;
            cur_parent = parent_edge->SrcNum();
            if(cur_parent == v_src)
                return true;
        }
        while(!tree_.IsRoot(cur_parent));
        return false;
    }

    bool ClonalGraphConstructor::VerticesBelongToSameTree(size_t v_src, size_t v_dst) const {
        return tree_.GetRootByVertex(v_dst) == tree_.GetRootByVertex(v_src);
    }

    bool ClonalGraphConstructor::ParentEdgeIsUndirected(size_t node) const {
        if(tree_.IsRoot(node))
            return false;
        auto parent_edge = tree_.GetParentEdge(node);
        return parent_edge->IsUndirected();
    }

    bool ClonalGraphConstructor::VerticesCanBeConnected(size_t v_src, size_t v_dst) const {
        if(tree_.IsIsolated(v_src) or tree_.IsIsolated(v_dst))
            return false;
        return v_src != v_dst and VerticesBelongToSameTree(v_src, v_dst) and
               !VerticesConnectedByDirectedPath(v_src, v_dst) and !ParentEdgeIsUndirected(v_src);
    }

    void ClonalGraphConstructor::AddDirectedEdges() {
        for(auto v1 = tree_.c_vertex_begin(); v1 != tree_.c_vertex_end(); v1++) {
            if(SkipVertex(*v1))
                continue;
            for (auto v2 = tree_.c_vertex_begin(); v2 != tree_.c_vertex_end(); v2++) {
                if(SkipVertex(*v2))
                    continue;
                if (clone_set_[*v1].VSHMs().size() >= clone_set_[*v2].VSHMs().size())
                    continue;
                // num SHMs in v1 < num SHMs in v2
                if (VerticesCanBeConnected(*v1, *v2)) {
                    //std::cout << "Good: " << *v1 << " " << *v2 << std::endl;
                    auto v_shms_1 = clone_set_[*v1].VSHMs();
                    auto v_shms_2 = clone_set_[*v2].VSHMs();
                    auto j_shms_1 = clone_set_[*v1].JSHMs();
                    auto j_shms_2 = clone_set_[*v2].JSHMs();
                    if (annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(v_shms_1, v_shms_2) and
                        annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(j_shms_1, j_shms_2)) {
                        clonal_graph_.AddNewEdge(*v1, *v2);
                    }
                }
            }
        }
    }

    //void ClonalGraphConstructor::InitializeClonalGraph() {
    //    AddEdgesFromTree();
    //    AddDirectedEdges();
    //    clonal_graph_.RemoveExtraEdges();
    //}
}