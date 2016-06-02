#pragma once

#include "cdr_labeling_primitives.hpp"
#include <read_archive.hpp>

namespace annotation_utils {
    enum StructuralRegion { CDR1, CDR2, CDR3 };

    std::ostream& operator<<(std::ostream& out, const StructuralRegion &region);

    class CDRAnnotatedClone {
        const core::Read &read_;

        std::unordered_map<StructuralRegion, seqan::Dna5String, std::hash<int>> region_string_map_;
        std::unordered_map<StructuralRegion, CDRRange, std::hash<int>> region_range_map_;

        void CheckRangeConsistencyFatal(CDRRange range);

        void UpdateStructuralRegion(StructuralRegion region, CDRRange range);

        void Initialize(CDRLabeling cdr_labeling);

    public:
        CDRAnnotatedClone(const core::Read &read, CDRLabeling cdr_labeling) : read_(read) {
            Initialize(cdr_labeling);
        }

        bool RegionExists(StructuralRegion region) const;

        CDRRange GetRangeByRegion(StructuralRegion region) const;

        seqan::Dna5String GetRegionString(StructuralRegion region) const;

        const core::Read& Read() const { return read_; }
    };

    std::ostream& operator<<(std::ostream& out, const CDRAnnotatedClone &obj);
}