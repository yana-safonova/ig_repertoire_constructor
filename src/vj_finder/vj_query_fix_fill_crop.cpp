#include <verify.hpp>
#include "vj_query_fix_fill_crop.hpp"

namespace vj_finder {
    VJHits AggressiveFillFixCropProcessor::Process(VJHits vj_hits) {
        core::Read read = vj_hits.Read();
        int right_shift = 0;
        int last_match_shift = 0;
        if(params_.crop_right and params_.fill_right) {
            auto j_hit = vj_hits.GetJHitByIndex(0);
            read.seq = seqan::prefix(read.seq, j_hit.LastMatchReadPos());
            auto j_suffix = seqan::suffix(j_hit.ImmuneGene().seq(), j_hit.LastMatchGenePos());
            read.seq += j_suffix;
            last_match_shift = int(seqan::length(j_suffix));
            if(params_.fix_right != 0) {
                for(size_t i = 0; i < params_.fix_right; i++) {
                    read.seq[read.length() - i - 1] = j_hit.ImmuneGene().seq()[j_hit.ImmuneGene().length() - i - 1];
                }
            }
        }
        int left_shift = 0;
        int first_match_shift = 0;
        if(params_.crop_left and params_.fill_right) {
            auto v_hit = vj_hits.GetVHitByIndex(0);
            read.seq = seqan::suffix(read.seq, v_hit.FirstMatchReadPos());
            seqan::Dna5String v_prefix = seqan::prefix(v_hit.ImmuneGene().seq(), v_hit.FirstMatchGenePos());
            v_prefix += read.seq;
            read.seq = v_prefix;
            left_shift = -int(v_hit.FirstMatchReadPos()) + int(v_hit.FirstMatchGenePos());
            first_match_shift = -int(v_hit.FirstMatchGenePos());
            if(params_.fix_left != 0) {
                for(size_t i = 0; i < params_.fix_left; i++) {
                    read.seq[i] = v_hit.ImmuneGene().seq()[i];
                }
            }
        }
        read_archive_.UpdateReadByIndex(read.id, read.seq);
        VJHits croppped_vj_hits(vj_hits.Read());
        croppped_vj_hits.AddVHit(vj_hits.GetVHitByIndex(0));
        croppped_vj_hits.AddJHit(vj_hits.GetJHitByIndex(0));
        croppped_vj_hits.AddLeftShift(left_shift);
        croppped_vj_hits.AddRightShift(left_shift + right_shift);
        croppped_vj_hits.ExtendFirstMatch(first_match_shift);
        croppped_vj_hits.ExtendLastMatch(last_match_shift);
        return croppped_vj_hits;
    }
}