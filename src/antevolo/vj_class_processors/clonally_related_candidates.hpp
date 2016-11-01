#pragma once

#include <annotation_utils/annotated_clone_set.hpp>

#include <unordered_map>

namespace antevolo {
    class ClonallyRelatedCandidates {
        const annotation_utils::CDRAnnotatedCloneSet &clone_set_;

        std::unordered_map<size_t, std::vector<size_t>> candidates_;
        std::set<size_t> all_elements_;

        void AddDirectedPair(size_t index1, size_t index2);

    public:
        ClonallyRelatedCandidates(const annotation_utils::CDRAnnotatedCloneSet &clone_set) :
                clone_set_(clone_set) { }

        void AddCandidatePair(size_t index1, size_t index2);

        std::vector<size_t> GetCandidatesForCloneIndex(size_t clone_index) const;

        bool Empty() const { return candidates_.size() == 0; }

        typedef std::set<size_t>::const_iterator CandidatesIterator;

        CandidatesIterator cbegin() const { return all_elements_.cbegin(); }

        CandidatesIterator cend() const { return all_elements_.cend(); }
    };
}