#pragma once

#include <memory>

#include "vj_alignment_structs.hpp"
#include "alignment_utils/alignment_positions.hpp"
#include "alignment_utils/pairwise_alignment.hpp"
#include "germline_utils/germline_gene_type.hpp"
#include "read_archive.hpp"
#include "vdj_alignments/hits_calculator/d_hits_calculator/abstract_d_gene_hits_calculator.hpp"
#include "d_gene_segment_hit.hpp"
#include "single_immune_gene_segment_hits.hpp"

namespace vdj_labeler {

class VDJHits {
private:
    const core::Read* read_ptr_;
    VGeneHits v_hits_;
    DGeneHits d_hits_;
    JGeneHits j_hits_;

public:
    VDJHits(const core::Read* read_ptr = nullptr);

    VDJHits(const core::Read* read_ptr,
            const std::vector<vj_finder::VGeneHit> &v_hits,
            const std::vector<vj_finder::JGeneHit> &j_hits);

    // VDJHits(const core::Read* read_ptr,
    //         const std::vector<vj_finder::VGeneHit> &v_hits,
    //         const std::vector<vj_finder::JGeneHit> &j_hits,
    //         AbstractDGeneHitsCalculator &d_gene_calculator);

    VDJHits(const vj_finder::VJHits &vj_hits);

    // VDJHits(const vj_finder::VJHits &vj_hits, AbstractDGeneHitsCalculator &d_gene_calculator);

    VDJHits(const VDJHits&)            = default;
    VDJHits(VDJHits&&)                 = default;
    VDJHits& operator=(const VDJHits&) = default;
    VDJHits& operator=(VDJHits&&)      = default;

    // void AddIgGeneAlignment(const germline_utils::SegmentType &gene_type,
    //                         alignment_utils::ImmuneGeneReadAlignment alignment);

    // void AddIgGeneAlignment(alignment_utils::ImmuneGeneReadAlignment alignment);

    size_t VHitsNumber() const { return v_hits_.size(); }
    size_t DHitsNumber() const { return d_hits_.size(); }
    size_t JHitsNumber() const { return j_hits_.size(); }

    const VGeneHits &VHits() const { return v_hits_; }
    const DGeneHits &DHits() const { return d_hits_; }
    const JGeneHits &JHits() const { return j_hits_; }

    void SetVHits(const VGeneHits &v_hits) { v_hits_ = v_hits; }
    void SetDHits(const DGeneHits &d_hits) { d_hits_ = d_hits; }
    void SetJHits(const JGeneHits &j_hits) { j_hits_ = j_hits; }

    const core::Read* ReadPtr() const { return read_ptr_; }
    const core::Read& Read() const {
        assert(read_ptr_ != nullptr);
        return *read_ptr_;
    }

    // const alignment_utils::ImmuneGeneReadAlignment& GetAlignmentByIndex(const germline_utils::SegmentType &gene_type,
    //                                                                     const size_t &index) const;
};

typedef std::shared_ptr<VDJHits> VDJHitsPtr;

} // End namespace vdj_labeler
