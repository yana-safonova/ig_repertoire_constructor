#pragma once

#include "vdj_alignments/vdj_hits.hpp"

namespace vdj_labeler {

class AbstractDGeneHitsCalculator {
protected:
    const core::ReadArchive &read_archive_;
    const germline_utils::ImmuneGeneDatabase &d_gene_database_;
    vdj_labeler::AlignmentQualityChecker &quality_checker_;

public:
    AbstractDGeneHitsCalculator(const core::ReadArchive &read_archive,
                                const germline_utils::ImmuneGeneDatabase &d_gene_database,
                                vdj_labeler::AlignmentQualityChecker &quality_checker) :
        read_archive_(read_archive),
        d_gene_database_(d_gene_database),
        quality_checker_(quality_checker)
    { }

    virtual ImmuneGeneSegmentHitsPtr ComputeDHits(core::ReadPtr read_ptr) = 0;

    virtual ~AbstractDGeneHitsCalculator() { }
};

} // End namespace vdj_labeler