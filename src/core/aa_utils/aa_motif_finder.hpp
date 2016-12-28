#pragma once

#include "amino_acid_motif.hpp"

#include <seqan/sequence.h>

namespace core {
    class AminoAcidMotifFinder {
    protected:
        AminoAcidMotifs aa_motifs_;

    public:
        AminoAcidMotifFinder(AminoAcidMotifs aa_motifs) :
                aa_motifs_(aa_motifs) { }

        virtual size_t FindMotif(seqan::String<seqan::AminoAcid> aa_seq, size_t pos);

        ~AminoAcidMotifFinder() { }
    };
}