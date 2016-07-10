#include "evolutionary_graph_constructor.hpp"

namespace antevolo {
    EvolutionaryEdgeConstructor* EvolutionaryGraphConstructor::GetEdgeConstructor() {
        return new SimpleEvolutionaryEdgeConstructor(config_.edge_construction_params);
    }
    /*
    EvolutionaryGraph EvolutionaryGraphConstructor::ConstructGraph() {
        EvolutionaryGraph graph;
        auto edge_constructor = GetEdgeConstructor();
        for(auto it = candidates_.cbegin(); it != candidates_.cend(); it++) {
            auto cur_candidates = candidates_.GetCandidatesForCloneIndex(*it);
            for(auto it2 = cur_candidates.begin(); it2 != cur_candidates.end(); it2++) {
                EvolutionaryEdge edge = edge_constructor->ConstructEdge(clone_set_[*it], clone_set_[*it2]);
                if(!edge.Empty())
                    graph.AddEdge(edge);
            }
        }
        TRACE("Graph contains " << graph.NumEdges() << " edges");
        TRACE("# " << EvolutionaryEdgeType::DirectedEdgeType << "s: " <<
                     graph.NumEdgesOfType(EvolutionaryEdgeType::DirectedEdgeType));
        TRACE("# " << EvolutionaryEdgeType::UndirectedEdgeType << "s: " <<
             graph.NumEdgesOfType(EvolutionaryEdgeType::UndirectedEdgeType));
        TRACE("# " << EvolutionaryEdgeType::IntersectedEdgeType << "s: " <<
             graph.NumEdgesOfType(EvolutionaryEdgeType::IntersectedEdgeType));
        //for(auto it = graph.cbegin(); it != graph.cend(); it++)
        //    std::cout << *it << std::endl;
        return graph;
    }
    */
}