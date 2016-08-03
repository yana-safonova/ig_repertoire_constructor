#include <verify.hpp>
#include "vj_query_fix_fill_crop.hpp"

namespace vj_finder {
    VJHits AggressiveFillFixCropProcessor::Process(VJHits vj_hits) {
        core::Read read = vj_hits.Read();
        int right_shift = 0;
        if(params_.crop_right and params_.fill_right) {
            auto j_hit = vj_hits.GetJHitByIndex(0);
            read.seq = seqan::prefix(read.seq, j_hit.LastMatchReadPos());
            auto j_suffix = seqan::suffix(j_hit.ImmuneGene().seq(), j_hit.LastMatchGenePos());
            read.seq += j_suffix;
            // todo: shift left bound of the first match
        }
        int left_shift = 0;
        if(params_.crop_left and params_.fill_right) {
            auto v_hit = vj_hits.GetVHitByIndex(0);
            read.seq = seqan::suffix(read.seq, v_hit.FirstMatchReadPos());
            seqan::Dna5String v_prefix = seqan::prefix(v_hit.ImmuneGene().seq(), v_hit.FirstMatchGenePos());
            v_prefix += read.seq;
            read.seq = v_prefix;
            left_shift = -int(v_hit.FirstMatchReadPos()) + int(v_hit.FirstMatchGenePos());
            // todo: shift right bound of the last match
        }
        read_archive_.UpdateReadByIndex(read.id, read.seq);
        vj_hits.AddLeftShift(left_shift);
        vj_hits.AddRightShift(left_shift + right_shift);
        return vj_hits;
    }
}