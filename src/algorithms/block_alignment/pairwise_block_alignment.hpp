#pragma once

#include <seqan/seq_io.h>
#include <seqan/align.h>

#include "block_alignment_primitives.hpp"
#include "block_alignment_utils.hpp"

namespace algorithms {

    struct PairwiseBlockAlignment {
        int kp_coverage;
        int int_score;
        AlignmentPath path;
        int start, finish;
        size_t needle_index;
        size_t needle_length;
        int overlap_length;
        double score;

        size_t first_match_read_pos() const {
            return path.first().read_pos;
        }

        size_t first_match_needle_pos() const {
            return path.first().needle_pos;
        }

        size_t last_match_needle_pos() const {
            return path.last().needle_pos + path.last().length;
        }

        size_t last_match_read_pos() const {
            return path.last().read_pos + path.last().length;
        }

        size_t left_half_segment_length() const {
            return last_match_needle_pos();
        }

        size_t right_half_segment_length() const {
            return needle_length - first_match_needle_pos();
        }

        size_t segment_length() const {
            return last_match_needle_pos() - first_match_needle_pos();
        }

        TAlign seqan_alignment(const seqan::Dna5String &read,
                               const seqan::Dna5String &gene) const {
            return path2seqanAlignment(this->path, read, gene);
        }

        std::string visualize(const seqan::Dna5String &read,
                              const seqan::Dna5String &gene) const {
            return this->path.visualize_matches(seqan::length(gene), seqan::length(read));
        }

        static PairwiseBlockAlignment path2Alignment(AlignmentPath &path,
                                        const seqan::Dna5String &read,
                                        const seqan::Dna5String &query,
                                        size_t needle_index,
                                        int score);
    };
}