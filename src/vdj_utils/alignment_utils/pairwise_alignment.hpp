#pragma once

#include "alignment_positions.hpp"

#include <seqan/align.h>

namespace alignment_utils {
    // Versatile pairwise alignment
    // T1 is a type of subject (e.g., ImmuneGene)
    // T2 is a type of query (e.g., Read)
    template<typename SubjectTypename, typename QueryTypename>
    class PairwiseAlignment {
        const SubjectTypename &subject_;
        const QueryTypename &query_;
        //AlignmentPositions positions_; it is important to have alignment positions?
        seqan::Align<seqan::Dna5String, seqan::ArrayGaps> alignment_;

        // computed characteristics
        size_t alignment_length_;
        size_t subject_alignment_length_;
        size_t query_alignment_length_;

        size_t num_gaps_;
        size_t num_matches_;
        size_t num_mismatches_;

        void ComputeAlignmentLengths() {
            VERIFY_MSG(false, "Implement me!");
        }

        void ComputeAlignmentStats() {
            VERIFY_MSG(false, "Implement me!");
        }

    public:
        PairwiseAlignment(const SubjectTypename &subject,
                          const QueryTypename &query,
                          seqan::Align<seqan::Dna5String, seqan::ArrayGaps> alignment) :
                subject_(subject), query_(query), alignment_(alignment) {
            //ComputeAlignmentLengths();
            //ComputeAlignmentStats();
        }

        const SubjectTypename & subject() const { return subject_; }

        const QueryTypename & query() const { return query_; }


        typedef seqan::Align<seqan::Dna5String, seqan::ArrayGaps> DnaGappedAlignment;

        const DnaGappedAlignment& Alignment() const { return alignment_; }

        size_t AlignmentLength() const { return alignment_length_; }

        size_t SubjectAlignmentLength() const { return subject_alignment_length_; }

        size_t QueryAlignmentLength() const { return query_alignment_length_; }

        bool Empty() const { return AlignmentLength() == 0; }


        size_t NumberGaps() const { return num_gaps_; }

        size_t NumberMatches() const { return num_matches_; }

        size_t NumberMismatches() const { return num_mismatches_; }


        size_t QueryPositionBySubjectPosition(size_t subject_pos) const {
            auto subject_row = seqan::row(alignment_, 0);
            auto query_row = seqan::row(alignment_, 1);
            size_t alignment_pos = seqan::toViewPosition(subject_row, subject_pos);
            return seqan::toSourcePosition(query_row, alignment_pos);
        }
    };

    typedef PairwiseAlignment<germline_utils::ImmuneGene, core::Read> ImmuneGeneReadAlignment;
}