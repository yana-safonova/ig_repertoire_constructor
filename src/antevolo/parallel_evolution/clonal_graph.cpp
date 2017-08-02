#include "clonal_graph.hpp"

namespace antevolo {
    void ClonalGraph::AddOldEdge(size_t src, size_t dst) {
        old_edges_[dst] = src;
        vertices_.insert(src);
        vertices_.insert(dst);
    }

    void ClonalGraph::AddNewEdge(size_t src, size_t dst) {
        if(new_edges_.find(dst) == new_edges_.end())
            new_edges_[dst] = std::set<size_t>();
        new_edges_[dst].insert(src);
    }

    /*
    void ClonalGraph::AddEdgesFromRhombSide(const ParallelRhomb &rhomb, RhombSide rhomb_side) {
        for (auto e = rhomb.edge_begin(rhomb_side); e != rhomb.edge_end(rhomb_side); e++) {
            size_t src = (*e)->SrcNum();
            size_t dst = (*e)->DstNum();
            AddEdge(src, dst);
        }
    }*/

    /*
    void ClonalGraph::InitializeEdges() {
        for(auto rhomb = stats_.rhombs.begin(); rhomb != stats_.rhombs.end(); rhomb++) {
            AddEdgesFromRhombSide(*rhomb, RhombSide1);
            AddEdgesFromRhombSide(*rhomb, RhombSide2);
        }
        for(auto it = tree_.cbegin(); it != tree_.cend(); it++) {
            auto edge_ptr = *it;
            if(!edge_ptr->IsDirected()) // ?
                continue;
            AddEdge(edge_ptr->SrcNum(), edge_ptr->DstNum());
        }
    }*/

    void ClonalGraph::FindTransitiveEdges() {
        for(auto it = new_edges_.begin(); it != new_edges_.end(); it++) {
            auto src_vertices = it->second;
            for(auto src1 = src_vertices.begin(); src1 != src_vertices.end(); src1++) {
                for(auto src2 = src_vertices.begin(); src2 != src_vertices.end(); src2++) {
                    if(*src1 == *src2)
                        continue;
                    if(old_edges_.find(*src1) != old_edges_.end()) {
                        if (old_edges_[*src1] == *src2) {
                            transitive_edges_.insert(std::make_pair(*src1, it->first));
                        }
                    }
                }
            }
        }
        //std::cout << "# transitive edges: " << transitive_edges_.size() << std::endl;
    }

    void ClonalGraph::RemoveTransitiveEdges() {
        for(auto e = transitive_edges_.begin(); e != transitive_edges_.end(); e++) {
            new_edges_[e->second].erase(e->first);
        }
        for(auto it = old_edges_.begin(); it != old_edges_.end(); it++) {
            if(all_edges_.find(it->second) == all_edges_.end())
                all_edges_[it->second] = std::set<size_t>();
            all_edges_[it->second].insert(it->first);
        }
        for(auto it = new_edges_.begin(); it != new_edges_.end(); it++) {
            auto src_vertices = it->second;
            for(auto src = src_vertices.begin(); src != src_vertices.end(); src++) {
                if(all_edges_.find(*src) == all_edges_.end())
                    all_edges_[*src] = std::set<size_t>();
                all_edges_[*src].insert(it->first);
            }
        }
    }

    void ClonalGraph::FindEndOfBulges() {
        std::map<size_t, size_t> vertex_mult;
        for(auto it = all_edges_.begin(); it != all_edges_.end(); it++) {
            auto outgoing_v = it->second;
            for(auto v = outgoing_v.begin(); v != outgoing_v.end(); v++) {
                if(vertex_mult.find(*v) == vertex_mult.end())
                    vertex_mult[*v] = 0;
                vertex_mult[*v]++;
            }
        }
        for(auto it = all_edges_.begin(); it != all_edges_.end(); it++) {
            auto outgoing_v = it->second;
            for (auto v = outgoing_v.begin(); v != outgoing_v.end(); v++) {
                if(vertex_mult[*v] > 1)
                    end_of_bulges_.insert(std::make_pair(it->first, *v));
            }
        }
    }

    void ClonalGraph::FillConflictingEdges() {
        for(auto it = all_edges_.begin(); it != all_edges_.end(); it++) {
            size_t src = it->first;
            auto dst_set = it->second;
            for(auto v = dst_set.begin(); v != dst_set.end(); v++) {
                if(!EdgeIsEndOfBulge(src, *v))
                    continue;
                if(conflicting_edges_.find(*v) == conflicting_edges_.end())
                    conflicting_edges_[*v] = std::set<size_t>();
                conflicting_edges_[*v].insert(src);
            }
        }
    }

    void ClonalGraph::RemoveExtraEdges() {
        FindTransitiveEdges();
        RemoveTransitiveEdges();
        FindEndOfBulges();
        FillConflictingEdges();
    }

    bool ClonalGraph::EdgeIsEndOfBulge(size_t src, size_t dst) const {
        return end_of_bulges_.find(std::make_pair(src, dst)) != end_of_bulges_.end();
    }

    std::set<size_t> ClonalGraph::OutgoingVertices(size_t src) const {
        VERIFY(all_edges_.find(src) != all_edges_.end());
        return all_edges_.at(src);
    }
}