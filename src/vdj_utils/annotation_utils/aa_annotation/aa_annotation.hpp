#pragma once

#include <verify.hpp>
#include <read_archive.hpp>
#include <seqan/translation.h>

namespace annotation_utils {
    typedef seqan::String<seqan::AminoAcid> AAString;

    template<class SequenceTypename>
    class AminoAcidAnnotation {
        const SequenceTypename* seq_prt_;
        AAString aa_seq_;
        size_t orf_;
        bool has_stop_codon_;
        bool in_frame_;

    public:
        AminoAcidAnnotation(const SequenceTypename& seq,
                            AAString aa_seq,
                            size_t orf,
                            bool has_stop_codon,
                            bool in_frame) : seq_prt_(&seq),
                                             aa_seq_(aa_seq),
                                             orf_(orf),
                                             has_stop_codon_(has_stop_codon),
                                             in_frame_(in_frame) { }

        const SequenceTypename& Seq() const { return *seq_prt_; }

        AAString AA() const { return aa_seq_; }

        size_t ORF() const { return orf_; }

        bool HasStopCodon() const { return has_stop_codon_; }

        bool InFrame() const { return in_frame_; }

        char GetAminoAcidByPos(size_t nucl_pos) const {
            VERIFY_MSG(nucl_pos < seq_prt_->length(), "Nucleotide position exceeds read length");
            if(nucl_pos < orf_)
                return '-';
            size_t aa_pos = (nucl_pos - orf_) / 3;
            if (aa_pos < seqan::length(aa_seq_))
                return aa_seq_[aa_pos];
            else if (aa_pos == seqan::length(aa_seq_))
                return '-';
            VERIFY_MSG(false, "Unknown position " << aa_pos << " of amino acid sequence " << aa_seq_);
            return '-';
        }
    };
}