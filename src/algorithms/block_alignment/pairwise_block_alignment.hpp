#pragma once

#include <seqan/seq_io.h>
#include <seqan/align.h>

#include "block_alignment_primitives.hpp"
#include "block_alignment_utils.hpp"

namespace algorithms {
    struct PairwiseBlockAlignment {
    private:
        int start_;
        int finish_;

    public:
        size_t subject_length;
        size_t query_length;

        size_t kp_coverage;
        int int_score;
        AlignmentPath path;

        int overlap_length;
        double score;

        int read_shift;

        PairwiseBlockAlignment() :
                start_(),
                finish_(),
                subject_length(),
                query_length(),
                kp_coverage(),
                int_score(),
                path(),
                overlap_length(),
                score(),
                read_shift() { }

        PairwiseBlockAlignment(AlignmentPath &path,
                               size_t subject_length,
                               size_t query_length,
                               int score);

        size_t first_match_read_pos() const {
            return size_t(path.first().read_pos);
        }

        size_t first_match_subject_pos() const {
            return size_t(path.first().subject_pos);
        }

        size_t last_match_subject_pos() const {
            return path.last().subject_pos + path.last().length;
        }

        size_t last_match_read_pos() const {
            return path.last().read_pos + path.last().length;
        }

        size_t left_half_segment_length() const {
            return last_match_subject_pos();
        }

        size_t right_half_segment_length() const {
            return subject_length - first_match_subject_pos();
        }

        size_t segment_length() const {
            return last_match_subject_pos() - first_match_subject_pos();
        }

        int start() const { return start_ + read_shift; }

        int finish() const { return finish_ + read_shift; }

        std::string visualize() const {
            return this->path.visualize_matches(static_cast<int>(subject_length), static_cast<int>(query_length));
        }

        void add_read_shift(int read_shift) {
            path.add_read_shift(read_shift);
            this->read_shift += read_shift;
        }

        void extend_first_match(int left_shift) {
            path.extend_first_match(left_shift);
        }

        void extend_last_match(int right_shift) {
            path.extend_last_match(right_shift);
        }
    };

    template<typename SubjectDatabase>
    class BlockAlignmentHits {
        const SubjectDatabase& db_;

        std::vector<std::pair<PairwiseBlockAlignment, size_t>> alignment_indices_;
        bool sorted_;

    public:
        BlockAlignmentHits(const SubjectDatabase& db) :
                db_(db), sorted_(false) { }

        BlockAlignmentHits(BlockAlignmentHits&& obj) : db_(obj.db_) {
            alignment_indices_ = std::move(obj.alignment_indices_);
            sorted_ = std::move(obj.sorted_);
        }

        BlockAlignmentHits& operator=(const BlockAlignmentHits& obj) {
            VERIFY(&db_ == &obj.db_);
            alignment_indices_ = obj.alignment_indices_;
            sorted_ = obj.sorted_;
            return *this;
        }

        void Add(PairwiseBlockAlignment block_alignment, size_t db_index) {
            VERIFY(db_index < db_.size()); // use Helper?
            alignment_indices_.push_back(std::make_pair(block_alignment, db_index));
        }

        typedef std::vector<std::pair<PairwiseBlockAlignment, size_t>>::const_iterator
                IndicedPairwiseBlockAlignmentConstIter;

        IndicedPairwiseBlockAlignmentConstIter cbegin() const { return alignment_indices_.cbegin(); }

        IndicedPairwiseBlockAlignmentConstIter cend() const { return alignment_indices_.cend(); }

        typedef std::vector<std::pair<PairwiseBlockAlignment, size_t>>::iterator
                IndicedPairwiseBlockAlignmentIter;

        IndicedPairwiseBlockAlignmentIter begin() { return alignment_indices_.begin(); }

        IndicedPairwiseBlockAlignmentIter end() { return alignment_indices_.end(); }

        size_t size() const { return alignment_indices_.size(); }

        typedef std::pair<PairwiseBlockAlignment, size_t> IndicedPairwiseBlockAlignment;

        const IndicedPairwiseBlockAlignment& operator[](size_t index) const {
            VERIFY(index < size());
            return alignment_indices_[index];
        }

        void SelectTopRecords(size_t limit = size_t(-1)) {
            limit = std::min(limit, size());
            if (limit == 0)
                return;
            using ctuple_type = decltype(*alignment_indices_.cbegin());
            auto score_function = [](const ctuple_type &a) { return a.first.int_score; };
            auto key_function = [&score_function](const ctuple_type &a) { return std::make_pair(-static_cast<int>(score_function(a)), a.second); };
            auto comp = [&key_function](const ctuple_type &a,
                                        const ctuple_type &b) -> bool { return key_function(a) < key_function(b); };
            // Return top <limit> positions
            std::nth_element(alignment_indices_.begin(), alignment_indices_.begin() + limit, alignment_indices_.end(), comp);
            alignment_indices_.resize(std::min(alignment_indices_.size(), limit));
            std::sort(alignment_indices_.begin(), alignment_indices_.end(), comp);
            sorted_ = true;
        }

        int BestScore() {
            if(!sorted_)
                SelectTopRecords();
            if(size() == 0)
                return 0;
            return alignment_indices_[0].first.int_score;
        }
    };
}