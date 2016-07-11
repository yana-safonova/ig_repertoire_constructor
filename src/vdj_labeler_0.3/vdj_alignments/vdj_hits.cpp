#include <memory>

#include "vdj_hits.hpp"
#include "immune_gene_alignment_converter.hpp"

namespace vdj_labeler {

ImmuneGeneSegmentHits::ImmuneGeneSegmentHits(const germline_utils::SegmentType &segment_type,
                                             const core::ReadPtr &read_ptr,
                                             const std::vector<vj_finder::ImmuneGeneHitPtr>& hits) :
    ImmuneGeneSegmentHits(segment_type, read_ptr)
{
    vj_finder::ImmuneGeneAlignmentConverter converter;
    for (auto& hit : hits) {
        assert(read_ptr != nullptr);
        auto conv_hit = converter.ConvertToAlignment(hit->ImmuneGene(), hit->Read(), hit->BlockAlignment(), true);
        hits_.push_back(std::make_shared<decltype(conv_hit)>(conv_hit));
    }
}

void ImmuneGeneSegmentHits::AddHit(const alignment_utils::ImmuneGeneReadAlignmentPtr &hit) {
    assert(hit->subject().GeneType().Segment() == segment_type_);
    hits_.push_back(hit);
}

alignment_utils::ImmuneGeneReadAlignmentPtr ImmuneGeneSegmentHits::operator[](const size_t &index) {
    assert(index < size());
    return hits_[index];
}

//-------------------------------------------------------------------------------

VDJHits::VDJHits(const core::ReadPtr &read_ptr,
                 const std::vector<vj_finder::ImmuneGeneHitPtr>& v_hits,
                 const std::vector<vj_finder::ImmuneGeneHitPtr>& j_hits) :
    read_ptr_(read_ptr),
    v_hits_(germline_utils::SegmentType::VariableSegment, read_ptr, v_hits),
    d_hits_(germline_utils::SegmentType::DiversitySegment, read_ptr),
    j_hits_(germline_utils::SegmentType::JoinSegment, read_ptr, j_hits)
{ }

void VDJHits::AddIgGeneAlignment(const alignment_utils::ImmuneGeneReadAlignmentPtr &alignment_ptr) {
    germline_utils::SegmentType segment_type = alignment_ptr->subject().GeneType().Segment();
    AddIgGeneAlignment(segment_type, alignment_ptr);
}

void VDJHits::AddIgGeneAlignment(const germline_utils::SegmentType &segment_type,
                                 const alignment_utils::ImmuneGeneReadAlignmentPtr &alignment_ptr) {
    if (segment_type == germline_utils::SegmentType::VariableSegment)
        v_hits_.AddHit(alignment_ptr);
    else if (segment_type == germline_utils::SegmentType::DiversitySegment)
        d_hits_.AddHit(alignment_ptr);
    else if (segment_type == germline_utils::SegmentType::JoinSegment)
        j_hits_.AddHit(alignment_ptr);
}

alignment_utils::ImmuneGeneReadAlignmentPtr VDJHits::GetAlignmentByIndex(
        const germline_utils::SegmentType &segment_type,
        const size_t &index)
{
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