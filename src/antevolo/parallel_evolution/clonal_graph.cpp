#include "clonal_graph.hpp"

namespace antevolo {
    void ClonalGraph::AddEdge(size_t src, size_t dst) {
        if (all_edges_.find(src) == all_edges_.end())
            all_edges_[src] = std::set<size_t>();
        all_edges_[src].insert(dst);
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
        for (auto src = all_edges_.begin(); src != all_edges_.end(); src++) { // a -> {b, c, d}
            auto outgoing = src->second; // {b, c, d}
            for (auto dst = outgoing.begin(); dst != outgoing.end(); dst++) { // b, ...
                if (all_edges_.find(*dst) == all_edges_.end())
                    continue;
                auto outgoing2 = all_edges_[*dst]; // b -> {c, e, f}
                for (auto dst2 = outgoing2.begin(); dst2 != outgoing2.end(); dst2++) { // c, ...
                    if (outgoing.find(*dst2) != outgoing.end()) // is c in {b ,c ,d}?
                        transitive_edges_.insert(std::make_pair(src->first, *dst2));
                }
            }
        }
        //std::cout << "# transitive edges: " << transitive_edges_.size() << std::endl;
    }

    void ClonalGraph::RemoveTransitiveEdges() {
        for(auto e = transitive_edges_.begin(); e != transitive_edges_.end(); e++) {
            if(all_edges_.find(e->first) != all_edges_.end()) {
                all_edges_[e->first].erase(e->second);
            }
        }
        //size_t num_edges = 0;
        //for(auto it = all_edges_.begin(); it != all_edges_.end(); it++)
        //    num_edges += it->second.size();
        //std::cout << "# edges after removing transitive ones: " << num_edges << std::endl;
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