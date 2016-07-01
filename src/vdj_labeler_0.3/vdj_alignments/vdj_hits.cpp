#include "vdj_hits.hpp"

namespace vdj_labeler {

void ImmuneGeneSegmentHits::AddHit(alignment_utils::ImmuneGeneReadAlignmentPtr hit) {
    assert(hit->subject().GeneType().Segment() == segment_type_);
    hits_.push_back(hit);
}

alignment_utils::ImmuneGeneReadAlignmentPtr ImmuneGeneSegmentHits::operator[](size_t index) {
    assert(index < size());
    return hits_[index];
}

//-------------------------------------------------------------------------------
void VDJHits::AddIgGeneAlignment(alignment_utils::ImmuneGeneReadAlignmentPtr alignment_ptr) {
    germline_utils::SegmentType segment_type = alignment_ptr->subject().GeneType().Segment();
    AddIgGeneAlignment(segment_type, alignment_ptr);
}

void VDJHits::AddIgGeneAlignment(germline_utils::SegmentType segment_type,
                                 alignment_utils::ImmuneGeneReadAlignmentPtr alignment_ptr) {
    if (segment_type == germline_utils::SegmentType::VariableSegment)
        v_hits_.AddHit(alignment_ptr);
    else if (segment_type == germline_utils::SegmentType::DiversitySegment)
        d_hits_.AddHit(alignment_ptr);
    else if (segment_type == germline_utils::SegmentType::JoinSegment)
        j_hits_.AddHit(alignment_ptr);
}

alignment_utils::ImmuneGeneReadAlignmentPtr VDJHits::GetAlignmentByIndex(germline_utils::SegmentType segment_type,
                                                                         size_t index) {
    if (segment_type == germline_utils::SegmentType::VariableSegment) {
        assert(index < v_hits_.size());
        return v_hits_[index];
    }
    else if (segment_type == germline_utils::SegmentType::DiversitySegment) {
        assert(index < d_hits_.size());
        return d_hits_[index];
    }
    assert(index < j_hits_.size());
    return j_hits_[index];
}

} // End namespace vdj_labeler