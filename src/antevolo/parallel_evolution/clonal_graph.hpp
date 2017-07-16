#pragma once

//#include "parallel_evolution_finder.hpp"

#include "../clone_set_with_fakes.hpp"
#include "../evolutionary_graph_utils/evolutionary_tree.hpp"

namespace antevolo {
    class ClonalGraph {
        //const EvolutionaryTree &tree_;
        //const CloneSetWithFakes &clone_set_;
        //const ParallelEvolutionStats &stats_;

        std::map<size_t, std::set<size_t> > all_edges_; // key - src vertex, values - dst vertices
        std::set<std::pair<size_t, size_t> > transitive_edges_;

        std::set<std::pair<size_t, size_t> > end_of_bulges_; // vertices with two or more incoming edges
        std::map<size_t, std::set<size_t> > conflicting_edges_; // key - dst vertex, values - src vertices

        //void AddEdgesFromRhombSide(const ParallelRhomb &rhomb, RhombSide rhomb_side);

        //void InitializeEdges();

        void FindTransitiveEdges();

        void RemoveTransitiveEdges();

        void FindEndOfBulges();

        void FillConflictingEdges();

    public:
        ClonalGraph(/*const EvolutionaryTree &tree ,
                    const ParallelEvolutionStats &stats*/) //: tree_(tree),
                                                           //clone_set_(tree.GetCloneSet())//,
                                                           //stats_(stats)
        {
            //InitializeEdges();
            //FindTransitiveEdges();
            //RemoveTransitiveEdges();
            //FindEndOfBulges();
            //FillConflictingEdges();
        }

        void AddEdge(size_t src, size_t dst);

        void RemoveExtraEdges();

        typedef std::map<size_t, std::set<size_t> >::const_iterator ClonalEdgeConstIterator;

        ClonalEdgeConstIterator cbegin() const { return all_edges_.cbegin(); }

        ClonalEdgeConstIterator cend() const { return all_edges_.cend(); }


        ClonalEdgeConstIterator conf_edge_cbegin() const { return conflicting_edges_.cbegin(); }

        ClonalEdgeConstIterator conf_edge_cend() const { return conflicting_edges_.cend(); }

        std::set<size_t> GetConflictingEdges(size_t dst) const {
            VERIFY(conflicting_edges_.find(dst) != conflicting_edges_.end());
            return conflicting_edges_.at(dst);
        }


        bool EdgeIsEndOfBulge(size_t src, size_t dst) const;

        std::set<size_t> OutgoingVertices(size_t src) const;

        //const CloneSetWithFakes& CloneSet() const { return clone_set_; }
    };
}