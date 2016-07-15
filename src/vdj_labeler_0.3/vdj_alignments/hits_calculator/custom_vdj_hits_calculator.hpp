#pragma once

#include "vj_finder/vj_alignment_info.hpp"
#include "vdj_alignments/aligners/gene_segment_aligner.hpp"

#include "vdj_hits_calculator.hpp"
#include "abstract_d_gene_hits_calculator.hpp"

namespace vdj_labeler {

class CustomVDJHitsCalculator: public VDJHitsCalculator {
    const VJAlignmentInfo &vj_alignment_info_;
    IgGeneHitsCalculator &v_hits_calculator_;
    IgGeneHitsCalculator &d_hits_calculator_;
    IgGeneHitsCalculator &j_hits_calculator_;

    void AddHits(VDJHitsPtr vdj_hits, IgGeneSegmentHitsPtr ig_gene_hits);

public:
    CustomVDJHitsCalculator(const FastqReadArchive &read_archive,
                            const VJAlignmentInfo &vj_alignment_info,
                            IgGeneHitsCalculator &v_hits_calculator,
                            IgGeneHitsCalculator &d_hits_calculator,
                            IgGeneHitsCalculator &j_hits_calculator) :
        VDJHitsCalculator(read_archive),
        vj_alignment_info_(vj_alignment_info),
        v_hits_calculator_(v_hits_calculator),
        d_hits_calculator_(d_hits_calculator),
        j_hits_calculator_(j_hits_calculator) { }

    VDJHitsStoragePtr ComputeHits();
};

} // End namespace vdj_labeler