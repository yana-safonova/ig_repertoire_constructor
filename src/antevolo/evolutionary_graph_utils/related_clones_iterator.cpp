#include "related_clones_iterator.hpp"
namespace antevolo {


    bool ClonesSharingCDR3sIterator::HasNext() {
        if (state_ == State::similarCDR3 && similar_cdr3s_it_ == similar_cdr3s_it_end_) {
            state_ = State::isDone;
        }
        return state_ != State::isDone;
    }

    const std::vector<size_t>& ClonesSharingCDR3sIterator::Next() {
        if (state_ == State::sameCDR3) {
            const auto& res = hamming_graph_info_.GetClonesByOldIndex(old_index_);
            state_ = State::similarCDR3;
            return res;
        }
        size_t neighbour_old_index = hamming_graph_info_.GetOldIndexByNewIndex(*similar_cdr3s_it_);
        const auto& res = hamming_graph_info_.GetClonesByOldIndex(neighbour_old_index);
        similar_cdr3s_it_++;
        return res;
    }

    bool RelatedClonesIterator::HasNext() {
        return current_clone_it_ != current_clone_it_end_ || vectors_iterator_.HasNext();
    }
    size_t RelatedClonesIterator::Next() {
        if (current_clone_it_ != current_clone_it_end_) {
            size_t next = *current_clone_it_;
            current_clone_it_++;
            return next;
        }
        const std::vector<size_t>& clones_sharing_current_cdr3 = vectors_iterator_.Next();
        current_clone_it_ = clones_sharing_current_cdr3.cbegin();
        current_clone_it_end_ = clones_sharing_current_cdr3.cend();
        size_t next = *current_clone_it_;
        current_clone_it_++;
        return next;
    }

    RelatedClonesIterator getRelatedClonesIterator(CDR3HammingGraphComponentInfo& hamming_graph_info, const std::string& cdr3) {
        size_t old_index = hamming_graph_info.GetOldIndexByCDR3(cdr3);
        ClonesSharingCDR3sIterator vectors_iterator(hamming_graph_info,
                                                    old_index);
        return RelatedClonesIterator(vectors_iterator);
    }

    RelatedClonesIterator getRelatedClonesIterator(CDR3HammingGraphComponentInfo& hamming_graph_info,
                                                   const annotation_utils::AnnotatedClone& clone) {
        std::string cdr3;
        const auto& cdr3_dna5 = clone.CDR3();
        size_t cdr3_length = seqan::length(cdr3_dna5);
        for (size_t i = 0; i < cdr3_length; ++i) {
            cdr3.push_back(cdr3_dna5[i]);
        }
        return getRelatedClonesIterator(hamming_graph_info, cdr3);
    }

}