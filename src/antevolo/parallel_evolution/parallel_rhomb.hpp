#pragma once

#include "../evolutionary_graph_utils/evolutionary_tree.hpp"

namespace antevolo {
    enum RhombSide { RhombSide1, RhombSide2 };

    class ParallelRhomb {
        const annotation_utils::CDRAnnotatedCloneSet& clone_set_;
        const EvolutionaryTree &tree_;

        std::vector<EvolutionaryEdgePtr> side1_edges_;
        std::vector<EvolutionaryEdgePtr> side2_edges_;

        std::vector<std::vector<annotation_utils::SHM> > side1_shms_;
        std::vector<std::vector<annotation_utils::SHM> > side2_shms_;


        std::vector<EvolutionaryEdgePtr>& GetEdgesBySide(RhombSide side);

        std::vector<std::vector<annotation_utils::SHM>>& GetSHMsBySide(RhombSide side);

    public:
        ParallelRhomb(const annotation_utils::CDRAnnotatedCloneSet& clone_set,
                      const EvolutionaryTree &tree) : clone_set_(clone_set),
                                                      tree_(tree) { }

        void AddEdgeToFrontOfSide(RhombSide side, EvolutionaryEdgePtr edge);

        void AddEdgeToEndOfSide(RhombSide side, EvolutionaryEdgePtr edge);


        size_t Side1Length() const { return side1_edges_.size(); }

        size_t Side2Length() const { return side2_edges_.size(); }

        size_t SideLength(RhombSide side) const { return (side == RhombSide1) ? Side1Length() : Side2Length(); }


        typedef std::vector<EvolutionaryEdgePtr>::const_iterator EvolutionaryEdgeIterator;

        EvolutionaryEdgeIterator edge_begin(RhombSide side) const;

        EvolutionaryEdgeIterator edge_end(RhombSide side) const;

        const EvolutionaryEdgePtr& GetEdgeByIndex(RhombSide side, size_t index) const;


        typedef std::vector<std::vector<annotation_utils::SHM> >::const_iterator SHMsIterator;

        SHMsIterator shms_begin(RhombSide side) const;

        SHMsIterator shms_end(RhombSide side) const;

        const std::vector<annotation_utils::SHM>& GetSHMsByIndex(RhombSide side, size_t index) const;


        size_t GetNumParallelSHMsBySide(RhombSide side) const;

        bool Ambiguous() const;

        size_t MinimalNumberParallelSHMs() const;

        //RhombSide GetParallelSide() const;


        const annotation_utils::CDRAnnotatedCloneSet& CloneSet() const { return clone_set_; }
    };

    std::ostream& operator<<(std::ostream &out, const ParallelRhomb &rhomb);
}
