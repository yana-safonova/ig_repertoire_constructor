#include "parallel_rhomb.hpp"

namespace antevolo {
    std::vector<EvolutionaryEdgePtr>& ParallelRhomb::GetEdgesBySide(RhombSide side) {
        if(side == RhombSide::RhombSide1)
            return side1_edges_;
        return side2_edges_;
    }

    std::vector<std::vector<annotation_utils::SHM> >& ParallelRhomb::GetSHMsBySide(RhombSide side) {
        if(side == RhombSide::RhombSide1)
            return side1_shms_;
        return side2_shms_;
    }

    void ParallelRhomb::AddEdgeToEndOfSide(RhombSide side, EvolutionaryEdgePtr edge) {
        auto &edges = GetEdgesBySide(side);
        edges.push_back(edge);
        auto &shms = GetSHMsBySide(side);
        auto added_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set_[edge->SrcNum()].VSHMs(),
                                                                        clone_set_[edge->DstNum()].VSHMs());
        shms.push_back(added_shms);
    }

    void ParallelRhomb::AddEdgeToFrontOfSide(RhombSide side, EvolutionaryEdgePtr edge) {
        auto &edges = GetEdgesBySide(side);
        edges.insert(edges.begin(), edge);
        auto &shms = GetSHMsBySide(side);
        auto added_shms = annotation_utils::SHMComparator::GetAddedSHMs(clone_set_[edge->SrcNum()].VSHMs(),
                                                                        clone_set_[edge->DstNum()].VSHMs());
        shms.insert(shms.begin(), added_shms);
    }

    ParallelRhomb::EvolutionaryEdgeIterator ParallelRhomb::edge_begin(RhombSide side) const {
        if(side == RhombSide::RhombSide1)
            return side1_edges_.begin();
        return side2_edges_.begin();
    }

    ParallelRhomb::EvolutionaryEdgeIterator ParallelRhomb::edge_end(RhombSide side) const {
        if(side == RhombSide::RhombSide1)
            return side1_edges_.end();
        return side2_edges_.end();
    }

    ParallelRhomb::SHMsIterator ParallelRhomb::shms_begin(RhombSide side) const {
        if(side == RhombSide::RhombSide1)
            return side1_shms_.begin();
        return side2_shms_.begin();
    }

    ParallelRhomb::SHMsIterator ParallelRhomb::shms_end(RhombSide side) const {
        if(side == RhombSide::RhombSide1)
            return side1_shms_.end();
        return side2_shms_.end();
    }

    const EvolutionaryEdgePtr& ParallelRhomb::GetEdgeByIndex(RhombSide side, size_t index) const {
        if(side == RhombSide::RhombSide1) {
            VERIFY(index < Side1Length());
            return side1_edges_[index];
        }
        VERIFY(index < Side2Length());
        return side2_edges_[index];
    }

    const std::vector<annotation_utils::SHM>& ParallelRhomb::GetSHMsByIndex(RhombSide side, size_t index) const {
        if(side == RhombSide::RhombSide1) {
            VERIFY(index < Side1Length());
            return side1_shms_[index];
        }
        VERIFY(index < Side2Length());
        return side2_shms_[index];
    }

    size_t ParallelRhomb::GetNumParallelSHMsBySide(RhombSide side) const {
        size_t parallel_shms = 0;
        if(side == RhombSide::RhombSide1) {
            for(size_t i = 0; i < Side1Length() - 1; i++)
                parallel_shms += side1_shms_[i].size();
            return parallel_shms;
        }
        for(size_t i = 0; i < Side2Length() - 1; i++)
            parallel_shms += side2_shms_[i].size();
        return parallel_shms;
    }

    bool ParallelRhomb::Ambiguous() const {
        return GetNumParallelSHMsBySide(RhombSide::RhombSide1) ==
                GetNumParallelSHMsBySide(RhombSide::RhombSide2);
    }

    size_t ParallelRhomb::MinimalNumberParallelSHMs() const {
        return std::min<size_t>(GetNumParallelSHMsBySide(RhombSide::RhombSide1),
                                GetNumParallelSHMsBySide(RhombSide::RhombSide2));
    }

    void print_edge_and_shms(std::ostream &out, const annotation_utils::CDRAnnotatedCloneSet &clone_set,
                             const EvolutionaryEdgePtr& edge,
                             const std::vector<annotation_utils::SHM>& shms) {
        out << "Edge " << edge->SrcNum() << "(" << clone_set[edge->SrcNum()].VSHMs().size() << ") -> " <<
            edge->DstNum() << "(" << clone_set[edge->DstNum()].VSHMs().size() << "): ";
        for(auto it = shms.cbegin(); it != shms.cend(); it++)
            out << *it << "; ";
    }

    std::ostream& operator<<(std::ostream &out, const ParallelRhomb &rhomb) {
        out << "Side 1: ";
        for(size_t i = 0; i < rhomb.Side1Length(); i++) {
            auto edge = rhomb.GetEdgeByIndex(RhombSide1, i);
            auto shms = rhomb.GetSHMsByIndex(RhombSide1, i);
            print_edge_and_shms(out, rhomb.CloneSet(), edge, shms);
        }
        out << std::endl;
        out << "Side 2: ";
        for(size_t i = 0; i < rhomb.Side2Length(); i++) {
            auto edge = rhomb.GetEdgeByIndex(RhombSide2, i);
            auto shms = rhomb.GetSHMsByIndex(RhombSide2, i);
            print_edge_and_shms(out, rhomb.CloneSet(), edge, shms);
        }
        out << std::endl;
        out << "Parallel SHMs: " << rhomb.GetNumParallelSHMsBySide(RhombSide1) << " vs " <<
            rhomb.GetNumParallelSHMsBySide(RhombSide2);
        return out;
    }
}
