#include "clonally_related_candidates.hpp"

namespace antevolo {
    void ClonallyRelatedCandidates::AddDirectedPair(size_t index1, size_t index2) {
        if(candidates_.find(index1) == candidates_.end())
            candidates_[index1] = std::vector<size_t>();
        candidates_[index1].push_back(index2);
    }

    void ClonallyRelatedCandidates::AddCandidatePair(size_t index1, size_t index2) {
        all_elements_.insert(index1);
        all_elements_.insert(index2);
        AddDirectedPair(index1, index2);
        AddDirectedPair(index2, index1);
    }

    std::vector<size_t> ClonallyRelatedCandidates::GetCandidatesForCloneIndex(size_t clone_index) const {
        if(candidates_.find(clone_index) == candidates_.end())
            return std::vector<size_t>();
        return candidates_.at(clone_index);
    }
}