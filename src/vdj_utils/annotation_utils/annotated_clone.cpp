#include <verify.hpp>

#include "annotated_clone.hpp"

#include <seqan/stream.h>
#include <seqan/translation.h>

namespace annotation_utils {
    std::ostream& operator<<(std::ostream& out, const StructuralRegion &region) {
        if(region == StructuralRegion::CDR1)
            out << "CDR1";
        else if(region == StructuralRegion::CDR2)
            out << "CDR2";
        else if(region == StructuralRegion::CDR3)
            out << "CDR3";
        else if(region == StructuralRegion::FR1)
            out << "FR1";
        else if(region == StructuralRegion::FR2)
            out << "FR2";
        else if(region == StructuralRegion::FR3)
            out << "FR3";
        else if(region == StructuralRegion::FR4)
            out << "FR4";
        else
            out << "Unknown region";
        return out;
    }

    void AnnotatedClone::CheckRangeConsistencyFatal(CDRRange range) {
        VERIFY(range.Full());
        VERIFY_MSG(range.start_pos < read_.length() and range.end_pos < read_.length(), "Start pos (" <<
                range.start_pos << ") or end pos (" << range.end_pos << ") exceeds read length " << read_.length());
    }

    void AnnotatedClone::UpdateStructuralRegion(StructuralRegion region, CDRRange range) {
        //TRACE("Updating " << region << " by range " << range);
        CheckRangeConsistencyFatal(range);
        region_range_map_[region] = range;
        seqan::Dna5String cdr_seq;
        if(range.Valid())
            cdr_seq = seqan::infixWithLength(read_.seq, range.start_pos, range.length());
        region_string_map_[region] = cdr_seq;
    }

    void AnnotatedClone::Initialize(CDRLabeling cdr_labeling) {
        UpdateStructuralRegion(StructuralRegion::FR1, CDRRange(0, cdr_labeling.cdr1.start_pos - 1));
        UpdateStructuralRegion(StructuralRegion::CDR1, cdr_labeling.cdr1);
        UpdateStructuralRegion(StructuralRegion::FR2, CDRRange(cdr_labeling.cdr1.end_pos + 1,
                                                               cdr_labeling.cdr2.start_pos - 1));
        UpdateStructuralRegion(StructuralRegion::CDR2, cdr_labeling.cdr2);
        UpdateStructuralRegion(StructuralRegion::FR3, CDRRange(cdr_labeling.cdr2.end_pos + 1,
                                                               cdr_labeling.cdr3.start_pos - 1));
        UpdateStructuralRegion(StructuralRegion::CDR3, cdr_labeling.cdr3);
        UpdateStructuralRegion(StructuralRegion::FR4, CDRRange(cdr_labeling.cdr3.end_pos + 1, read_.length() - 1));
    }

    bool AnnotatedClone::RegionIsEmpty(StructuralRegion region) const {
        if (region_range_map_.find(region) == region_range_map_.end())
            return true;
        return seqan::length(GetRegionString(region)) == 0;
    }

    seqan::Dna5String AnnotatedClone::GetRegionString(StructuralRegion region) const {
        if(region_range_map_.find(region) == region_range_map_.end())
            return seqan::Dna5String();
        //VERIFY_MSG(region_range_map_.find(region) != region_range_map_.end(),
        //           "Clone does not have information about region " << region);
        return region_string_map_.at(region);
    }

    CDRRange AnnotatedClone::GetRangeByRegion(StructuralRegion region) const {
        VERIFY_MSG(!RegionIsEmpty(region), "Clone does not have information about region " << region);
        return region_range_map_.at(region);
    }

    const alignment_utils::ImmuneGeneReadAlignment& AnnotatedClone::GetAlignmentBySegment(
            germline_utils::SegmentType segment_type) const {
        VERIFY_MSG(segment_type == germline_utils::SegmentType::VariableSegment or
                           segment_type == germline_utils::SegmentType::JoinSegment,
                   "Segment " << segment_type << "is not variable or join");
        if(segment_type == germline_utils::SegmentType::VariableSegment)
            return VAlignment();
        return JAlignment();
    }

    size_t AnnotatedClone::GetAminoAcidPosByNucleotidePos(size_t nucl_pos) const {
        VERIFY_MSG(nucl_pos < read_.length(), "Position " << nucl_pos << " exceeds sequence length");
        return (nucl_pos - ORF()) / 3;
    }

    char AnnotatedClone::GetAminoAcidByNucleotidePos(size_t nucl_pos) const {
        VERIFY_MSG(nucl_pos < read_.length(), "Position " << nucl_pos << " exceeds sequence length");
        return aa_annotation_.AA()[GetAminoAcidPosByNucleotidePos(nucl_pos)];
    }

    StructuralRegion AnnotatedClone::GetRegionBySHM(SHM shm) const {
        if(shm.read_nucl_pos < GetRangeByRegion(StructuralRegion::CDR1).start_pos)
            return StructuralRegion::FR1;
        if(shm.read_nucl_pos <= GetRangeByRegion(StructuralRegion::CDR1).end_pos)
            return StructuralRegion::CDR1;
        if(shm.read_nucl_pos < GetRangeByRegion(StructuralRegion::CDR2).start_pos)
            return StructuralRegion::FR2;
        if(shm.read_nucl_pos <= GetRangeByRegion(StructuralRegion::CDR2).end_pos)
            return StructuralRegion::CDR2;
        if(shm.read_nucl_pos < GetRangeByRegion(StructuralRegion::CDR3).start_pos)
            return StructuralRegion::FR3;
        if(shm.read_nucl_pos <= GetRangeByRegion(StructuralRegion::CDR3).end_pos)
            return StructuralRegion::CDR3;
        return StructuralRegion::FR4;
    }

    std::ostream& operator<<(std::ostream& out, const AnnotatedClone &obj) {
        out << obj.Read() << std::endl;
        if(!obj.RegionIsEmpty(StructuralRegion::CDR1))
            out << "CDR1: " << obj.GetRegionString(StructuralRegion::CDR1) << std::endl;
        if(!obj.RegionIsEmpty(StructuralRegion::CDR2))
            out << "CDR1: " << obj.GetRegionString(StructuralRegion::CDR2) << std::endl;
        if(!obj.RegionIsEmpty(StructuralRegion::CDR3))
            out << "CDR1: " << obj.GetRegionString(StructuralRegion::CDR3) << std::endl;
        return out;
    }
}