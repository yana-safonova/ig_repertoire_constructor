#pragma once

#include "alignment_utils/pairwise_alignment.hpp"

namespace recombination_utils {

class CleavedIgGeneAlignment {
    alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment_ptr_;
    int left_cleavage_length_;
    int right_cleavage_length_;
    int num_left_shms_;
    int num_right_shms_;
    size_t total_num_shms_;

public:
    CleavedIgGeneAlignment(alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment_ptr,
                           int left_cleavage_length,
                           int right_cleavage_length,
                           int num_left_shms,
                           int num_right_shms) :
            gene_alignment_ptr_(gene_alignment_ptr),
            left_cleavage_length_(left_cleavage_length),
            right_cleavage_length_(right_cleavage_length),
            num_left_shms_(num_left_shms),
            num_right_shms_(num_right_shms)
    {
        total_num_shms_ = size_t(int(gene_alignment_ptr_->NumberSHMs()) + num_left_shms_ + num_right_shms_);
    }

    CleavedIgGeneAlignment(const CleavedIgGeneAlignment& cleaved_gene) :
            gene_alignment_ptr_(cleaved_gene.gene_alignment_ptr_),
            left_cleavage_length_(cleaved_gene.left_cleavage_length_),
            right_cleavage_length_(cleaved_gene.right_cleavage_length_),
            num_left_shms_(cleaved_gene.num_left_shms_),
            num_right_shms_(cleaved_gene.num_right_shms_),
            total_num_shms_(cleaved_gene.total_num_shms_)
    { }

    alignment_utils::ImmuneGeneReadAlignmentPtr GeneAlignment() const { return gene_alignment_ptr_; }

    // negative length of cleavage means existence of the left palindrome of this length
    // positive length of cleavage shows length of gene cleavage
    int LeftCleavageLength() const { return left_cleavage_length_; }

    int RightCleavageLength() const { return right_cleavage_length_; }

    size_t SHMsNumber() const { return total_num_shms_; }

    // number of somatic mutations in left event
    int NumberLeftSHMs() const { return num_left_shms_; }

    // number of somatic mutations in right event
    int NumberRightSHMs() const { return num_right_shms_; }

    size_t GeneId() const { return gene_alignment_ptr_->subject().id(); }

    size_t StartReadPosition() const { return size_t(int(gene_alignment_ptr_->StartQueryPosition()) +
                                              left_cleavage_length_); }

    size_t EndReadPosition() const { return size_t(int(gene_alignment_ptr_->EndQueryPosition()) +
                                            right_cleavage_length_ * -1); }

    bool LeftEventIsEmpty() const { return left_cleavage_length_ == 0; }

    bool LeftEventIsPalindrome() const { return left_cleavage_length_ < 0; }

    bool LeftEventIsCleavage() const { return left_cleavage_length_ > 0; }

    bool RightEventIsEmpty() const { return right_cleavage_length_ == 0; }

    bool RightEventIsPalindrome() const { return right_cleavage_length_ < 0; }

    bool RightEventIsCleavage() const { return right_cleavage_length_ > 0; }
};

std::ostream& operator<<(std::ostream &out, const CleavedIgGeneAlignment& cleaved_gene);

} // End namespace recombination_utils