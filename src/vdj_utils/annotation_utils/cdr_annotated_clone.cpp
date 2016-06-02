#include "cdr_annotated_clone.hpp"

#include <seqan/stream.h>

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

    void CDRAnnotatedClone::CheckRangeConsistencyFatal(CDRRange range) {
        VERIFY(range.Full());
        VERIFY_MSG(range.start_pos < read_.length() and range.end_pos < read_.length(), "Start pos (" <<
                range.start_pos << ") or end pos (" << range.end_pos << ") exceeds read length " << read_.length());
    }

    void CDRAnnotatedClone::UpdateStructuralRegion(StructuralRegion region, CDRRange range) {
        CheckRangeConsistencyFatal(range);
        region_range_map_[region] = range;
        region_string_map_[region] = seqan::infixWithLength(read_.seq, range.start_pos, range.length());
    }

    void CDRAnnotatedClone::Initialize(CDRLabeling cdr_labeling) {
        UpdateStructuralRegion(StructuralRegion::CDR1, cdr_labeling.cdr1);
        UpdateStructuralRegion(StructuralRegion::CDR2, cdr_labeling.cdr2);
        UpdateStructuralRegion(StructuralRegion::CDR3, cdr_labeling.cdr3);
    }

    bool CDRAnnotatedClone::RegionExists(StructuralRegion region) const {
        return region_range_map_.find(region) != region_range_map_.end();
    }

    seqan::Dna5String CDRAnnotatedClone::GetRegionString(StructuralRegion region) const {
        VERIFY_MSG(RegionExists(region), "Clone does not have information about region " << region);
        return region_string_map_.at(region);
    }

    CDRRange CDRAnnotatedClone::GetRangeByRegion(StructuralRegion region) const {
        VERIFY_MSG(RegionExists(region), "Clone does not have information about region " << region);
        return region_range_map_.at(region);
    }

    std::ostream& operator<<(std::ostream& out, const CDRAnnotatedClone &obj) {
        out << obj.Read() << std::endl;
        if(obj.RegionExists(StructuralRegion::CDR1))
            out << "CDR1: " << obj.GetRegionString(StructuralRegion::CDR1) << std::endl;
        if(obj.RegionExists(StructuralRegion::CDR2))
            out << "CDR1: " << obj.GetRegionString(StructuralRegion::CDR2) << std::endl;
        if(obj.RegionExists(StructuralRegion::CDR3))
            out << "CDR1: " << obj.GetRegionString(StructuralRegion::CDR3) << std::endl;
        return out;
    }
}