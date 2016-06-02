#include "hcdr2_labeler.hpp"

#include "aa_utils/aa_motif_finder.hpp"

#include <seqan/translation.h>

namespace cdr_labeler {
    using namespace annotation_utils;

    size_t HCDR2Labeler::ComputeStartPosition(const germline_utils::ImmuneGene &immune_gene,
                                              CDRRange previous_cdr) {
        using namespace seqan;
        StringSet<String<AminoAcid>, Owner<ConcatDirect<> > > aa_seqs;
        translate(aa_seqs, immune_gene.seq(), SINGLE_FRAME);
        std::cout << "V gene aa seq: " << aa_seqs[0] << std::endl;
        core::AminoAcidMotifs aa_start_motif(params_.residues_before);
        //size_t aa_length = immune_gene.length() / 3;
        size_t start_range = previous_cdr.end_pos / 3 + 1 + params_.distance_from_cdr1_end -
                params_.distance_shift - aa_start_motif.length();
        size_t end_range = start_range + 2 * params_.distance_shift;
        size_t best_pos = size_t(-1);
        size_t best_score = size_t(-1);
        core::AminoAcidMotifFinder *aa_motif_finder = new core::AminoAcidMotifFinder(aa_start_motif);
        for(size_t i = start_range; i < end_range; i++) {
            size_t cur_score = aa_motif_finder->FindMotif(aa_seqs[0], i);
            if(cur_score < best_score) {
                best_pos = i * 3 + aa_start_motif.length() * 3;
                best_score = cur_score;
            }
        }
        return best_pos;
    }

    size_t HCDR2Labeler::ComputeEndPosition(const germline_utils::ImmuneGene &, size_t) {
        VERIFY_MSG(false, "Implement me!");
        return size_t(-1);
    }

    CDRRange HCDR2Labeler::ComputeRange(const germline_utils::ImmuneGene &immune_gene, CDRRange previous_cdr) {
        size_t start_pos = ComputeStartPosition(immune_gene, previous_cdr);
        size_t end_pos = ComputeEndPosition(immune_gene, start_pos);
        return CDRRange(start_pos, end_pos);
    }
}