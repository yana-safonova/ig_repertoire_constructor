#include <verify.hpp>

#include "aa_calculator.hpp"

namespace annotation_utils {
    bool SimpleAACalculator::ComputeInFrame(const CDRLabeling &cdr_labeling) const {
        CDRRange end_region = (cdr_labeling.cdr3.Valid()) ? cdr_labeling.cdr3 : cdr_labeling.cdr2;
        VERIFY_MSG(end_region.Valid() and cdr_labeling.cdr1.Valid(),
                   "CDRs regions are not defined, ORF cannot be identified");
        return (end_region.end_pos - cdr_labeling.cdr1.start_pos + 1) % 3 == 0;
    }

    bool SimpleAACalculator::FindStopCodon(const AAString &aa_str) const {
        bool has_stop_codon = false;
        for(size_t i = 0; i < seqan::length(aa_str); i++)
            if(aa_str[i] == '*') {
                has_stop_codon = true;
                break;
            }
        return has_stop_codon;
    }

    AminoAcidAnnotation<core::Read> SimpleAACalculator::ComputeAminoAcidAnnotation(const core::Read &read,
                                                                       const CDRLabeling &cdr_labeling) const {
        VERIFY_MSG(cdr_labeling.cdr1.Valid(), "CDR1 is not defined, AA sequence cannot be computed");
        using namespace seqan;
        StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > aa_seqs;
        size_t orf = cdr_labeling.cdr1.start_pos % 3;
        translate(aa_seqs, suffix(read.seq, orf), SINGLE_FRAME);
        VERIFY(seqan::length(aa_seqs) > 0);
        AAString aa_seq = aa_seqs[0];
        bool in_frame = ComputeInFrame(cdr_labeling);
        bool has_stop_codon = FindStopCodon(aa_seq);
        return AminoAcidAnnotation<core::Read>(read, aa_seq, orf, has_stop_codon, in_frame);
    }
}