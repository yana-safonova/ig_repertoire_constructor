#pragma once

#include "vj_alignment_info.hpp"
#include "vdj_alignments/aligners/gene_segment_aligner.hpp"
#include "vdj_hits_calculator.hpp"
#include "vdj_alignments/hits_calculator/d_hits_calculator/info_based_d_hits_calculator.hpp"
#include "vdj_alignments/vdj_hits.hpp"

namespace vdj_labeler {

class CustomVDJHitsCalculator: public VDJHitsCalculator {
    const vj_finder::VJAlignmentInfo &vj_alignment_info_;
    InfoBasedDHitsCalculator &d_hits_calculator_;

    void AddHits(VDJHitsPtr vdj_hits, ImmuneGeneSegmentHitsPtr ig_gene_hits);

public:
    CustomVDJHitsCalculator(const vj_finder::VJAlignmentInfo &vj_alignment_info,
                            InfoBasedDHitsCalculator &d_hits_calculator) :
        VDJHitsCalculator(),
        vj_alignment_info_(vj_alignment_info),
        d_hits_calculator_(d_hits_calculator)
    { }

    VDJHitsStoragePtr ComputeHits() const override;
};

} // End namespace vdj_labeler