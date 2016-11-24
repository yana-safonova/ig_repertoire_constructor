#include "related_clones_iterator.hpp"
namespace antevolo {


    bool ClonesSharingCDR3sIterator::HasNext() {
        if (state_ == State::similarCDR3 && similar_cdr3s_it_ == similar_cdr3s_it_end_) {
            state_ = State::isDone;
        }
        return state_ != State::isDone;
    }
    const std::vector<size_t>& ClonesSharingCDR3sIterator::operator*() const {
        if (state_ == State::sameCDR3) {
            return hamming_graph_info_.GetClonesByOldIndex(old_index_);
        }
        return hamming_graph_info_.GetClonesByOldIndex(*similar_cdr3s_it_);
    }
    ClonesSharingCDR3sIterator& ClonesSharingCDR3sIterator::operator++() {
        if (state_ == State::sameCDR3) {
            state_ = State::similarCDR3;
            return *this;
        }
        similar_cdr3s_it_++;
        return *this;
    }
    ClonesSharingCDR3sIterator ClonesSharingCDR3sIterator::operator++(int) {
        ClonesSharingCDR3sIterator res(*this);
        ++*this;
        return res;
    }

    bool RelatedClonesIterator::HasNext() {
        return current_vector_it_ != clones_sharing_current_cdr3_.cend() || vectors_iterator_.HasNext();
    }
    size_t RelatedClonesIterator::Next() {
        if (current_vector_it_ != clones_sharing_current_cdr3_.cend()) {
            size_t next = *current_vector_it_;
            current_vector_it_++;
            return next;
        }
        ++vectors_iterator_;
        current_vector_it_ = (*vectors_iterator_).cbegin();
        size_t next = *current_vector_it_;
        current_vector_it_++;
        return next;
    }

    RelatedClonesIterator getRelatedClonesIterator(CDR3HammingGraphInfo& hamming_graph_info, const std::string& cdr3) {
        size_t old_index = hamming_graph_info.GetOldIndexByCDR3(cdr3);
        ClonesSharingCDR3sIterator vectors_iterator(hamming_graph_info,
                                                    old_index);
        return RelatedClonesIterator(vectors_iterator);
    }

    RelatedClonesIterator getRelatedClonesIterator(CDR3HammingGraphInfo& hamming_graph_info,
                                                   const annotation_utils::AnnotatedClone& clone) {
        std::string cdr3;
        const auto& cdr3_dna5 = clone.CDR3();
        size_t  cdr3_length = seqan::length(cdr3);
        for (size_t i = 0; i < cdr3_length; ++i) {
            cdr3.push_back(cdr3_dna5[i]);
        }
        return getRelatedClonesIterator(hamming_graph_info, cdr3);
    }

}