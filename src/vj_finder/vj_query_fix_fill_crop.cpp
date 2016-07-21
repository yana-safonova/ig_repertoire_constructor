#include <verify.hpp>
#include "vj_query_fix_fill_crop.hpp"

namespace vj_finder {
    VJHits FillFixCropProcessor::PerformFixing(VJHits vj_hits) {
        if(!params_.enable_fixing)
            return vj_hits;
        core::Read read = vj_hits.Read();
        for(size_t i = 0; i < params_.fix_left; i++) {
            int gene_pos = static_cast<int>(i) - vj_hits.GetVHitByIndex(0).LeftShift();
            if(gene_pos >= 0 and gene_pos < static_cast<int>(vj_hits.GetVHitByIndex(0).ImmuneGene().length()))
                read.seq[i] = vj_hits.GetVHitByIndex(0).ImmuneGene().seq()[gene_pos];
        }
        for(size_t i = read.length() - params_.fix_right; i < read.length(); i++) {
            int gene_pos = static_cast<int>(i) - vj_hits.GetVHitByIndex(0).RightShift();
            if (gene_pos >= 0 && gene_pos < static_cast<int>(vj_hits.GetJHitByIndex(0).ImmuneGene().length()))
                read.seq[i] = vj_hits.GetJHitByIndex(0).ImmuneGene().seq()[gene_pos];
        }
        vj_hits.UpdateRead(read);
        return vj_hits;
//        VJHits fixed_vj_hits(read);
//        for(size_t i = 0; i < vj_hits.NumVHits(); i++) {
//            fixed_vj_hits.AddVHit(vj_hits.GetVHitByIndex(i));
//        }
//        for(size_t i = 0; i < vj_hits.NumJHits(); i++)
//            fixed_vj_hits.AddJHit(vj_hits.GetJHitByIndex(i));
//        return fixed_vj_hits;
    }

    VJHits FillFixCropProcessor::PerformFillingCropping(VJHits vj_hits) {
        TRACE("Read: " << vj_hits.Read().name << ", start: " << vj_hits.GetVHitByIndex(0).Start() << ", end: " <<
                     vj_hits.GetVHitByIndex(0).End());
        core::Read read = vj_hits.Read();
        TRACE("Right filling & cropping...");
        int right_shift = 0;
        if(params_.crop_right and vj_hits.GetJHitByIndex(0).End() < static_cast<int>(read.length())) {
            // no alignment editing in this case
            read.seq = seqan::prefix(read.seq, vj_hits.GetJHitByIndex(0).End());
        }
        else if(params_.fill_right and vj_hits.GetJHitByIndex(0).End() > static_cast<int>(read.length())) {
            auto j_gene = vj_hits.GetJHitByIndex(0).ImmuneGene();
            auto j_suffix = seqan::suffix(j_gene.seq(),
                                          j_gene.length() - (vj_hits.GetJHitByIndex(0).End() - read.length()));
            read.seq += j_suffix;
            // extend end alignment of J gene and read by length of j_suffix
            right_shift = int(seqan::length(j_suffix));
        }
        TRACE("Left filling & cropping...");
        int left_shift = 0;
        if(params_.crop_left and vj_hits.GetVHitByIndex(0).Start() > 0)  {
            // shift to left all alignment positions of read by vj_hits.GetVHitByIndex(0).Start()
            read.seq = seqan::suffix(read.seq, vj_hits.GetVHitByIndex(0).Start());
            left_shift = -vj_hits.GetVHitByIndex(0).Start();
        }
        else if(params_.fill_left and vj_hits.GetVHitByIndex(0).Start() < 0) {
            seqan::Dna5String germline_prefix = seqan::prefix(vj_hits.GetVHitByIndex(0).ImmuneGene().seq(),
                                                              -vj_hits.GetVHitByIndex(0).Start());
            left_shift = int(seqan::length(germline_prefix));
            germline_prefix += read.seq;
            read.seq = germline_prefix;
            // shift to right all alignment positions of read by length of germline_prefix
        }
        vj_hits.UpdateRead(read);
        vj_hits.AddLeftShift(left_shift);
        vj_hits.AddRightShift(left_shift + right_shift);
        //VJHits filled_cropped_vj_hits(read);
        //VGeneHit v_hit = vj_hits.GetVHitByIndex(0);
        //v_hit.AddShift(left_shift);
        //filled_cropped_vj_hits.AddVHit(v_hit);
        //JGeneHit j_hit = vj_hits.GetJHitByIndex(0);
        //j_hit.AddShift(left_shift + right_shift);
        //filled_cropped_vj_hits.AddJHit(j_hit);
        TRACE("Start: " << vj_hits.GetVHitByIndex(0).Start() << ", end: " <<
                      vj_hits.GetVHitByIndex(0).End());
        return vj_hits;
    }

    VJHits FillFixCropProcessor::Process(VJHits vj_hits) {
        TRACE("Fixing, filling and cropping read " << vj_hits.Read() << " starts");
        VJHits fixed_vj_hits = PerformFixing(vj_hits);
        TRACE("Fixing done");
        VJHits final_vj_hits = PerformFillingCropping(fixed_vj_hits);
        TRACE("Filling and croppping done");
        return final_vj_hits;
    }
}