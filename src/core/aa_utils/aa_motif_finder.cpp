#include "aa_motif_finder.hpp"

namespace core {
    // computation hamming distance between motif and aa sequence
    size_t AminoAcidMotifFinder::FindMotif(seqan::String<seqan::AminoAcid> aa_seq, size_t pos) {
        size_t best_score = 0;
        for(auto it = aa_motifs_.cbegin(); it != aa_motifs_.cend(); it++) {
            auto aa_motif = *it;
            if(pos + aa_motif.length() >= seqan::length(aa_seq))
                continue;
            size_t cur_score = 0;
            for(size_t i = 0; i < aa_motif.length(); i++) {
                if(aa_motif[i] == char(aa_seq[pos + i]) and aa_motif[i] != 'X')
                    cur_score++;
            }
            if(cur_score > best_score)
                best_score = cur_score;
        }
        return best_score;
    }
}