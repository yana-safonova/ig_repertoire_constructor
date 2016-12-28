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
        else
            out << "Unknown region";
        return out;
    }

    void AnnotatedClone::CheckRangeConsistencyFatal(CDRRange range) {
        VERIFY(range.Full());
        VERIFY_MSG(range.start_pos < read_ptr_->length() and range.end_pos < read_ptr_->length(),
                   "Start pos (" << range.start_pos << ") or end pos (" <<
                           range.end_pos << ") exceeds read length " << read_ptr_->length());
    }

    void AnnotatedClone::UpdateStructuralRegion(StructuralRegion region, CDRRange range) {
        //TRACE("Updating " << region << " by range " << range);
        CheckRangeConsistencyFatal(range);
        region_range_map_[region] = range;
        seqan::Dna5String cdr_seq;
        if(range.Valid())
            cdr_seq = seqan::infixWithLength(read_ptr_->seq, range.start_pos, range.length());
        region_string_map_[region] = cdr_seq;
    }

    void AnnotatedClone::Initialize(CDRLabeling cdr_labeling) {
        UpdateStructuralRegion(StructuralRegion::CDR1, cdr_labeling.cdr1);
        UpdateStructuralRegion(StructuralRegion::CDR2, cdr_labeling.cdr2);
        UpdateStructuralRegion(StructuralRegion::CDR3, cdr_labeling.cdr3);
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