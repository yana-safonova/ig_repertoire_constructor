#include <verify.hpp>
#include "annotated_evolutionary_tree.hpp"

namespace antevolo {
    void AnnotatedEvolutionaryTree::CheckConsistencyFatal(size_t clone_id) const {
        VERIFY_MSG(tree_.ContainsClone(clone_id), "Evolutionary tree does not contain clone " << clone_id);
    }

    void AnnotatedEvolutionaryTree::InitializeCloneSHMMap() {
        for(auto it = shm_map_.c_shm_clone_begin(); it != shm_map_.c_shm_clone_end(); it++) {
            auto edges = it->second;
            for(auto edge = edges.begin(); edge != edges.end(); edge++) {
                if(clone_added_shm_map_.find(edge->second) == clone_added_shm_map_.end())
                    clone_added_shm_map_[edge->second] = 0;
                clone_added_shm_map_[edge->second]++;
            }
        }
        size_t root_id = tree_.GetRoot();
        clone_added_from_root_map_[root_id] = 0;
        std::queue<size_t> child_queue;
        child_queue.push(root_id);
        while(!child_queue.empty()) {
            size_t cur_vertex = child_queue.front();
            child_queue.pop();
            if(tree_.IsLeaf(cur_vertex))
                continue;
            auto outgoing_edges = tree_.OutgoingEdges(cur_vertex);
            for(auto e = outgoing_edges.begin(); e != outgoing_edges.end(); e++) {
                size_t dst = (*e)->DstNum();
                clone_added_from_root_map_[dst] = clone_added_shm_map_[dst] + clone_added_from_root_map_[cur_vertex];
                child_queue.push(dst);
            }
        }
    }

    size_t AnnotatedEvolutionaryTree::RootDepth() const {
        size_t root_id = tree_.GetRoot();
        return clone_set_[root_id].VSHMs().size() + clone_set_[root_id].JSHMs().size();
    }

    size_t AnnotatedEvolutionaryTree::TreeDepth() const {
        size_t tree_depth = 0;
        for(auto it = tree_.c_vertex_begin(); it != tree_.c_vertex_end(); it++) {
            if(tree_.IsLeaf(*it)) {
                size_t num_added_shms = clone_added_from_root_map_.at(*it);
                tree_depth = std::max(tree_depth, num_added_shms);
            }
        }
        return tree_depth;
    }

    size_t AnnotatedEvolutionaryTree::GetRegionLength(annotation_utils::StructuralRegion region) const {
        size_t root_id = tree_.GetRoot();
        return seqan::length(clone_set_[root_id].GetRegionString(region));
    }
}