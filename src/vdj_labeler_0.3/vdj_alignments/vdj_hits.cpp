#include <memory>

#include "vdj_hits.hpp"
#include "immune_gene_alignment_converter.hpp"

namespace vdj_labeler {

void ImmuneGeneSegmentHits::AddHit(alignment_utils::ImmuneGeneReadAlignmentPtr hit) {
    assert(hit->subject().GeneType().Segment() == segment_type_);
    hits_.push_back(hit);
}

alignment_utils::ImmuneGeneReadAlignmentPtr ImmuneGeneSegmentHits::operator[](size_t index) {
    assert(index < size());
    return hits_[index];
}

ImmuneGeneSegmentHits::ImmuneGeneSegmentHits(germline_utils::SegmentType segment_type, core::ReadPtr read_ptr,
                                             const std::vector<vj_finder::ImmuneGeneHitPtr>& hits) :
    ImmuneGeneSegmentHits(segment_type, read_ptr)
{
    vj_finder::ImmuneGeneAlignmentConverter converter;
    for (auto& hit : hits) {
        assert(read_ptr != nullptr);
        // TODO change *read_ptr to hit->Read(). Now hit->Read() supposingly has a bug.
        auto conv_hit = converter.ConvertToAlignment(hit->ImmuneGene(), *read_ptr, hit->BlockAlignment());
        hits_.push_back(std::make_shared<decltype(conv_hit)>(conv_hit));
    }
}

//-------------------------------------------------------------------------------

VDJHits::VDJHits(const core::ReadPtr read_ptr,
                 const std::vector<vj_finder::ImmuneGeneHitPtr>& v_hits,
                 const std::vector<vj_finder::ImmuneGeneHitPtr>& j_hits) :
    read_ptr_(read_ptr),
    v_hits_(germline_utils::SegmentType::VariableSegment, read_ptr, v_hits),
    d_hits_(germline_utils::SegmentType::DiversitySegment, read_ptr),
    j_hits_(germline_utils::SegmentType::JoinSegment, read_ptr, j_hits)
{ }

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