#include "logger/logger.hpp"
#include "right_event_shms_calculator.hpp"

using namespace vdj_labeler;

int RightEventSHMsCalculator::ComputeNumberCleavedSHMs(
        const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
        const size_t cleavage_length) const
{
    TRACE("Computation of # SHMs in right cleavage of length " << cleavage_length);
    // We have a test for this function
    size_t alignment_cleavage = gene_alignment->SubjectAlignmentLength() - 1 - gene_alignment->EndSubjectPosition();
    assert(cleavage_length >= alignment_cleavage);
    size_t rel_cleavage_length = cleavage_length - alignment_cleavage;

    auto alignment = gene_alignment->Alignment();
    auto &gene = seqan::row(alignment, 0);
    auto &read = seqan::row(alignment, 1);
    int num_shms = 0;
    int cur_cleavage = 0;
    for (int i = static_cast<int>(seqan::length(gene)) - 1; i >= 0; i--) {
        if (read[i] != '-')
            cur_cleavage++;
        if (gene[i] != read[i])
            num_shms++;
        if (cur_cleavage == int(rel_cleavage_length))
            break;
    }
    TRACE("Cleavage length: " << cleavage_length << ", rel cleavage length: " << rel_cleavage_length);
    TRACE("#SHMs: -" << num_shms);
    return -1 * num_shms;
}

int RightEventSHMsCalculator::ComputeNumberPalindromeSHMs(
        const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
        const size_t palindrome_length) const
{
    TRACE("Computation of #SHMs in right palindrome of length " << palindrome_length);
    // if gene has alignment to read with gaps at the end, we can not compute palindrome
    assert(gene_alignment->EndSubjectPosition() == gene_alignment->GeneLength() - 1);
    int num_shms = 0;
    for (size_t i = 0; i < palindrome_length; i++) {
        size_t gene_pos = gene_alignment->SubjectAlignmentLength() - 1 - i;
        size_t read_pos = gene_alignment->EndQueryPosition() + 1 + i;
        seqan::DnaString nucl(gene_alignment->subject().seq()[gene_pos]);
        seqan::complement(nucl);
        // nucl has len 1.
        if (nucl[0] != gene_alignment->query().seq[read_pos])
            num_shms++;
    }
    TRACE("#SHMs: +" << num_shms);
    return num_shms;
}

int RightEventSHMsCalculator::ComputeNumberSHMs(
        const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
        const int,
        const int right_cleavage_length) const
{
    if (right_cleavage_length == 0)
        return 0;
    // if gene was cleaved, number of SHMs would 0 or negative
    // since some unaligned nucleotides were cleaved
    if (right_cleavage_length > 0)
        return ComputeNumberCleavedSHMs(gene_alignment, size_t(right_cleavage_length));
    // if gene contains palindrome, number of SHMs would be 0 or positive
    // since some nucleotides in palindrome are mutated
    return ComputeNumberPalindromeSHMs(gene_alignment, size_t(right_cleavage_length * -1));
}

int RightEventSHMsCalculator::ComputeNumberSHMsForRightEvent(
        const alignment_utils::ImmuneGeneReadAlignmentPtr gene_alignment,
        const int right_cleavage_length) const
{
    if (right_cleavage_length == 0)
        return 0;
    // if gene was cleaved, number of SHMs would 0 or negative
    // since some unaligned nucleotides were cleaved
    if (right_cleavage_length > 0)
        return ComputeNumberCleavedSHMs(gene_alignment, size_t(right_cleavage_length));
    // if gene contains palindrome, number of SHMs would be 0 or positive
    // since some nucleotides in palindrome are mutated
    return ComputeNumberPalindromeSHMs(gene_alignment, size_t(right_cleavage_length * -1));
}