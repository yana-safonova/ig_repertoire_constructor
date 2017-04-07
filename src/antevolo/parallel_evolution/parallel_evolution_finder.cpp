#include "parallel_evolution_finder.hpp"

namespace antevolo {
    bool ParallelEvolutionFinder::VerticesConnectedByDirectedPath(size_t v_src, size_t v_dst) const {
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

    bool ParallelEvolutionFinder::VerticesBelongToSameTree(size_t v_src, size_t v_dst) const {
        return tree_.GetRootByVertex(v_dst) == tree_.GetRootByVertex(v_src);
    }

    bool ParallelEvolutionFinder::ParentEdgeIsUndirected(size_t node) const {
        if(tree_.IsRoot(node))
            return false;
        auto parent_edge = tree_.GetParentEdge(node);
        return parent_edge->IsUndirected();
    }

    bool ParallelEvolutionFinder::PairIsGood(size_t v_src, size_t v_dst) const {
        if(tree_.IsIsolated(v_src) or tree_.IsIsolated(v_dst))
            return false;
        return v_src != v_dst and VerticesBelongToSameTree(v_src, v_dst) and
                !VerticesConnectedByDirectedPath(v_src, v_dst) and !ParentEdgeIsUndirected(v_src);
    }

    void ParallelEvolutionFinder::FillAddedDirectedEdges() {
        for(auto v1 = tree_.c_vertex_begin(); v1 != tree_.c_vertex_end(); v1++)
            for(auto v2 = tree_.c_vertex_begin(); v2 != tree_.c_vertex_end(); v2++) {
                if(clone_set_[*v1].VSHMs().size() >= clone_set_[*v2].VSHMs().size())
                    continue;
                // num SHMs in v1 < num SHMs in v2
                if(PairIsGood(*v1, *v2)) {
                    //std::cout << "Good: " << *v1 << " " << *v2 << std::endl;
                    auto v_shms_1 = clone_set_[*v1].VSHMs();
                    auto v_shms_2 = clone_set_[*v2].VSHMs();
                    if(annotation_utils::SHMComparator::SHMs1AreNestedInSHMs2(v_shms_1, v_shms_2)) {
                        if(added_directed_edges_map_.find(*v1) == added_directed_edges_map_.end()) {
                            added_directed_edges_map_[*v1] = std::set<size_t>();
                        }
                        added_directed_edges_map_[*v1].insert(*v2);
                    }
                }
            }
        PrintAddedNestedEdges();
    }

    void ParallelEvolutionFinder::RemoveUndirectedRhombs() {
        //INFO("Removing undirected rhombs");
        for(auto it = added_directed_edges_map_.begin(); it != added_directed_edges_map_.end(); it++) {
            std::vector<size_t> nodes_to_remove;
            for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
                if(tree_.IsRoot(*it2))
                    continue;
                auto parent_edge = tree_.GetParentEdge(*it2);
                if(parent_edge->IsUndirected()) // and it->second.find(parent_edge->SrcNum()) != it->second.end())
                    nodes_to_remove.push_back(*it2);
            }
            for(auto n = nodes_to_remove.begin(); n != nodes_to_remove.end(); n++)
                it->second.erase(*n);
        }
        //PrintAddedNestedEdges();
    }

    void ParallelEvolutionFinder::RemoveIncomingNestedRhombs() {
        //INFO("Removing incoming nested rhombs");
        std::vector<size_t> keys;
        std::vector<size_t> nodes_to_remove;
        for(auto it = added_directed_edges_map_.begin(); it != added_directed_edges_map_.end(); it++)
            keys.push_back(it->first);
        for(size_t i = 0; i < keys.size(); i++)
            for(size_t j = 0; j < keys.size(); j++) {
                if(i == j)
                    continue;
                size_t v1 = keys[i];
                size_t v2 = keys[j];
                if(VerticesConnectedByDirectedPath(v1, v2))
                    nodes_to_remove.push_back(v1);
            }
        //INFO("# nodes to be removed as parts of nested rhombs: " << nodes_to_remove.size());
        for(auto it = nodes_to_remove.begin(); it != nodes_to_remove.end(); it++)
            added_directed_edges_map_.erase(*it);
        //PrintAddedNestedEdges();
    }

    void ParallelEvolutionFinder::RemoveNestedOutgoingRhombs() {
        //INFO("Removing outgoing nested rhombs");
        size_t num_removed_nodes= 0;
        for(auto it = added_directed_edges_map_.begin(); it != added_directed_edges_map_.end(); it++) {
            std::vector<size_t> values;
            for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++)
                values.push_back(*it2);
            std::vector<size_t> nodes_to_remove;
            for(size_t i = 0; i < values.size(); i++)
                for(size_t j = 0; j < values.size(); j++) {
                    if(i == j)
                        continue;
                    size_t v1 = values[i];
                    size_t v2 = values[j];
                    if(VerticesConnectedByDirectedPath(v1, v2))
                        nodes_to_remove.push_back(v1);
                }
            num_removed_nodes += nodes_to_remove.size();
            for(auto n = nodes_to_remove.begin(); n != nodes_to_remove.end(); n++)
                it->second.erase(*n);
        }
        //PrintAddedNestedEdges();
    }

    void ParallelEvolutionFinder::PrintAddedNestedEdges() {
        size_t num_added_edges = 0;
        for(auto it = added_directed_edges_map_.begin(); it != added_directed_edges_map_.end(); it++)
            num_added_edges += it->second.size();
        //INFO("# nested edges: " << num_added_edges);
    }

    EvolutionaryEdgePtr ParallelEvolutionFinder::GetAddedEvolutionaryTree(size_t v_src, size_t v_dst) {
        return EvolutionaryEdgePtr(new UndirectedEvolutionaryEdge(clone_set_[v_src], clone_set_[v_dst], v_src, v_dst));
    }

    ParallelRhomb ParallelEvolutionFinder::ComputeRhombByAddedEdge(size_t v_src, size_t v_dst) {
        //INFO("Construction of rhomb on " << v_src << " -> " << v_dst);
        ParallelRhomb rhomb(clone_set_, tree_); // added edge v_src -> v_dst will be always a part of side1
        size_t common_source = v_src;
        rhomb.AddEdgeToFrontOfSide(RhombSide1, GetAddedEvolutionaryTree(v_src, v_dst));
        //INFO("Constructing side 1");
        while(!VerticesConnectedByDirectedPath(common_source, v_dst)) {
            auto parent_edge = tree_.GetParentEdge(common_source);
            //INFO("Adding edge " << parent_edge->SrcNum() << " -> " << parent_edge->DstNum());
            rhomb.AddEdgeToFrontOfSide(RhombSide1, parent_edge);
            common_source = parent_edge->SrcNum();
        }
        //INFO("Constructing side 2");
        size_t current_dst = v_dst;
        while(current_dst != common_source) {
            auto parent_edge = tree_.GetParentEdge(current_dst);
            //INFO("Adding edge " << parent_edge->SrcNum() << " -> " << parent_edge->DstNum());
            rhomb.AddEdgeToFrontOfSide(RhombSide2, parent_edge);
            current_dst = parent_edge->SrcNum();
        }
        return rhomb;
    }

    bool ParallelEvolutionFinder::RhombIsGood(ParallelRhomb rhomb) {
        std::vector<RhombSide> sides = {RhombSide1, RhombSide2};
        for(auto s = sides.cbegin(); s != sides.cend(); s++) {
            size_t num_non_zero = 0;
            for(size_t i = 0; i < rhomb.SideLength(*s); i++)
                if(rhomb.GetSHMsByIndex(*s, i).size() != 0)
                    num_non_zero++;
            if(num_non_zero == 1)
                return false;
        }
        return true;
    }

    std::vector<ParallelRhomb> ParallelEvolutionFinder::ComputeRhombs() {
        std::vector<ParallelRhomb> rhombs;
        for(auto it = added_directed_edges_map_.begin(); it != added_directed_edges_map_.end(); it++)
            for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
                size_t v_src = it->first;
                size_t v_dst = *it2;
                auto rhomb = ComputeRhombByAddedEdge(v_src, v_dst);
                if(!RhombIsGood(rhomb))
                    continue;
                std::cout << rhomb << std::endl;
                std::cout << "# parallel SHMs: " << rhomb.MinimalNumberParallelSHMs() << std::endl;
                std::cout << "------------" << std::endl;
                rhombs.push_back(rhomb);
            }
        return rhombs;
    }

    ParallelEvolutionStats ParallelEvolutionFinder::ComputeStats(const std::vector<ParallelRhomb> &rhombs) {
        ParallelEvolutionStats stats(rhombs.size());
        for(auto it = rhombs.begin(); it != rhombs.end(); it++)
            stats.AddParallelSHM(it->MinimalNumberParallelSHMs());
        return stats;
    }

    ParallelEvolutionStats ParallelEvolutionFinder::ComputeParallelSHMs() {
        //INFO("== Processing tree " << tree_.GetTreeOutputFname(""));
        FillAddedDirectedEdges();
        if(added_directed_edges_map_.size() == 0)
            return ParallelEvolutionStats();
        RemoveUndirectedRhombs();
        RemoveIncomingNestedRhombs();
        RemoveNestedOutgoingRhombs();
        //for(auto it = added_directed_edges_map_.begin(); it != added_directed_edges_map_.end(); it++)
        //    for(auto it2 = it->second.begin(); it2 != it->second.end(); it2++)
        //        std::cout << it->first << " - " << *it2 << ", # added SHMs: " <<
        //                (clone_set_[*it2].VSHMs().size() - clone_set_[it->first].VSHMs().size())<< std::endl;
        auto parallel_rhombs = ComputeRhombs();
        //INFO("# true parallel rhombs: " << parallel_rhombs.size());
        return ComputeStats(parallel_rhombs);
    }
}
