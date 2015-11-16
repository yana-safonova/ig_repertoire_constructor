#pragma once

#include "../vdj_alignments/alignment_structs.hpp"

class CleavedIgGeneAlignment {
    const IgGeneAlignment &gene_alignment_;
    int left_cleavage_length_;
    int right_cleavage_length_;

public:
    CleavedIgGeneAlignment(IgGeneAlignment gene_alignment,
                           int left_cleavage_length,
                           int right_cleavage_length) :
            gene_alignment_(gene_alignment),
            left_cleavage_length_(left_cleavage_length),
            right_cleavage_length_(right_cleavage_length) { }

    const IgGeneAlignment& GeneAlignment() const { return gene_alignment_; }

    // negative length of cleavage means existance of the left palindrome of this length
    // positive length of cleavage shows length of gene cleavage
    int LeftCleavageLength() const { return left_cleavage_length_; }

    int RightCleavageLength() const { return right_cleavage_length_; }
};